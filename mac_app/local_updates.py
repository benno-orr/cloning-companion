"""User-confirmed updates from builds published on this Mac (no remote feed)."""
from __future__ import annotations

import hashlib
import json
import os
import plistlib
import re
import subprocess
import sys
import tempfile
import time
import uuid
from pathlib import Path

BUNDLE_ID = "com.bennoorr.cloningcompanion"


def update_directory() -> Path:
    path = Path.home() / "Library" / "Application Support" / "CloningCompanion" / "Updates"
    path.mkdir(parents=True, exist_ok=True)
    return path


def atomic_json(path: Path, data) -> None:
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    temporary.write_text(json.dumps(data, indent=2), encoding="utf-8")
    temporary.chmod(0o600)
    temporary.replace(path)


def version_key(value: str) -> tuple:
    if not re.fullmatch(r"\d+\.\d+\.\d+", value):
        raise ValueError("Update version must be major.minor.patch")
    return tuple(map(int, value.split(".")))


def bundle_info(path: Path) -> dict:
    if path.is_symlink() or path.suffix != ".app" or not path.is_dir():
        raise ValueError("Update target must be a real .app bundle, not a symlink")
    with (path / "Contents/Info.plist").open("rb") as handle:
        info = plistlib.load(handle)
    if info.get("CFBundleIdentifier") != BUNDLE_ID:
        raise ValueError("This is not a CloningCompanion application")
    if info.get("CFBundleExecutable") != "CloningCompanion":
        raise ValueError("Unexpected application executable")
    version_key(info["CFBundleShortVersionString"])
    return info


def bundle_digest(path: Path) -> str:
    """Content checksum includes relative names, executable bits and symlinks."""
    digest = hashlib.sha256()
    for entry in sorted(path.rglob("*")):
        name = entry.relative_to(path).as_posix()
        if entry.is_symlink():
            # Framework links must stay within the signed bundle.
            entry.resolve().relative_to(path.resolve())
            digest.update(json.dumps([name, "link", os.readlink(entry)]).encode())
        elif entry.is_file():
            digest.update(json.dumps([name, "file", entry.stat().st_mode & 0o111, entry.stat().st_size]).encode())
            with entry.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
    return digest.hexdigest()


def verify_bundle(path: Path, manifest: dict) -> None:
    info = bundle_info(path)
    if info["CFBundleShortVersionString"] != manifest["version"]:
        raise ValueError("Build changed since publication; wait for the build to finish")
    if bundle_digest(path) != manifest["sha256"]:
        raise ValueError("Update checksum failed; the current application was not changed")
    subprocess.run(["/usr/bin/codesign", "--verify", "--deep", "--strict", str(path)],
                   check=True, capture_output=True, text=True)


def publish(path: Path, directory: Path) -> dict:
    path = path.resolve()
    info = bundle_info(path)
    manifest = {"version": info["CFBundleShortVersionString"], "bundle": str(path),
                "sha256": bundle_digest(path), "channel": "local-builds"}
    verify_bundle(path, manifest)
    directory.mkdir(parents=True, exist_ok=True)
    atomic_json(directory / "latest.json", manifest)
    return manifest


def running_bundle() -> Path | None:
    if not getattr(sys, "frozen", False):
        return None
    return next((p for p in Path(sys.executable).resolve().parents if p.suffix == ".app"), None)


def check_update(current: str, directory: Path, target: Path | None, manifest: dict | None = None) -> dict:
    manifest_path = directory / "latest.json"
    if manifest is None and not manifest_path.exists():
        return {"ok": True, "available": False, "message": "No local build has been published yet."}
    if manifest is None:
        manifest = json.loads(manifest_path.read_text())
    newer = version_key(manifest["version"]) > version_key(current)
    if not newer:
        return {"ok": True, "available": False, "message": f"You’re up to date with local builds ({current})."}
    if not Path(manifest["bundle"]).is_dir():
        raise ValueError("Published build is unavailable. Finish a new local build and check again.")
    reason = None
    if target is None:
        reason = "Open the desktop .app to install updates; source/browser previews cannot replace themselves."
    elif str(target).startswith("/Volumes/") or "AppTranslocation" in target.parts:
        reason = "Move the app to Applications before updating."
    return {"ok": True, "available": True, "version": manifest["version"],
            "canInstall": reason is None, "message": reason or "A new local build is ready to install."}


def stage_update(target: Path, directory: Path, current: str, session: dict, manifest: dict | None = None) -> Path:
    """Prepare and verify everything before the running app is asked to quit."""
    target = target.absolute()
    bundle_info(target)
    status = check_update(current, directory, target, manifest)
    if not status.get("available") or not status.get("canInstall"):
        raise ValueError(status["message"])
    if manifest is None:
        manifest = json.loads((directory / "latest.json").read_text())
    source = Path(manifest["bundle"])
    if source.resolve() == target.resolve():
        raise ValueError("The published build cannot replace itself")
    # A sibling staging directory guarantees same-filesystem atomic renames and
    # checks destination write access now, while the old app is still running.
    mode = "bundle"
    try:
        stage_root = Path(tempfile.mkdtemp(prefix=".CloningCompanion-update-", dir=target.parent))
    except PermissionError:
        # /Applications may be administrator-owned while this particular app
        # is user-owned. Replace its Contents without writing to /Applications.
        if not os.access(target, os.W_OK):
            raise PermissionError("This app is not writable. Ask an administrator to install the update.")
        mode = "contents"
        stage_root = Path(tempfile.mkdtemp(prefix="staged-", dir=directory))
        if stage_root.stat().st_dev != target.stat().st_dev:
            raise ValueError("App and update staging must share a filesystem for safe replacement")
    staged = stage_root / "CloningCompanion.app"
    subprocess.run(["/usr/bin/ditto", str(source), str(staged)], check=True, capture_output=True)
    verify_bundle(staged, manifest)
    if mode == "contents" and any({p.name for p in bundle.iterdir()} != {"Contents"} for bundle in (target, staged)):
        raise ValueError("This app layout requires an administrator-assisted replacement")
    job_dir = Path(tempfile.mkdtemp(prefix="install-", dir=directory))
    atomic_json(job_dir / "recovery.plasmidverify", {"format": "CloningCompanion Project", "version": 1, "state": session["state"]})
    atomic_json(job_dir / "session.json", session)
    backup_parent = target.parent if mode == "bundle" else stage_root
    config = {"target": str(target), "staged": str(staged), "parentPid": os.getpid(), "mode": mode,
              "backup": str(backup_parent / f"CloningCompanion.previous-{uuid.uuid4().hex[:8]}.app"),
              "manifest": manifest, "originalDigest": bundle_digest(target)}
    atomic_json(job_dir / "job.json", config)
    return job_dir / "job.json"


def replace_bundle(config: dict, launch=None) -> None:
    """Atomic replacement; restore the old bundle if replacement/launch fails."""
    target, staged, backup = (Path(config[key]) for key in ("target", "staged", "backup"))
    bundle_info(target)
    contents_only = config.get("mode") == "contents"
    allowed_parent = staged.parent if contents_only else target.parent
    if backup.exists() or backup.parent != allowed_parent or (not contents_only and staged.parent.parent != target.parent):
        raise ValueError("Unsafe update destinations")
    if target == staged or target == backup:
        raise ValueError("Update destinations overlap")
    if bundle_digest(target) != config["originalDigest"]:
        raise ValueError("Installed app changed while update was pending; nothing was replaced")
    verify_bundle(staged, config["manifest"])
    if contents_only:
        if staged.stat().st_dev != target.stat().st_dev or any(
                {p.name for p in bundle.iterdir()} != {"Contents"} for bundle in (target, staged)):
            raise ValueError("Cannot safely replace the app contents")
    if launch is None:
        launch = lambda path: subprocess.run(["/usr/bin/open", "-n", str(path)], check=True, capture_output=True)
    source_part, staged_part, backup_part = target, staged, backup
    if contents_only:
        backup.mkdir()
        source_part, staged_part, backup_part = target / "Contents", staged / "Contents", backup / "Contents"
    source_part.rename(backup_part)
    try:
        staged_part.rename(source_part)
        launch(target)
    except Exception:
        if source_part.exists():
            source_part.rename(staged_part)
        backup_part.rename(source_part)
        launch(target)
        raise


def apply_update_job(job_path: Path) -> int:
    """Runs in a separate instance of the staged frozen executable."""
    directory = job_path.parent
    config = json.loads(job_path.read_text())
    atomic_json(directory / "ready.json", {"ready": True})
    parent_exited = False
    try:
        deadline = time.monotonic() + 120
        while True:
            try:
                os.kill(int(config["parentPid"]), 0)
            except ProcessLookupError:
                parent_exited = True
                break
            if time.monotonic() > deadline:
                raise TimeoutError("App did not quit; the update was cancelled without replacing it")
            time.sleep(.2)
        # Recovery is available to the new app before it starts.
        atomic_json(update_directory() / "pending-session.json", {"job": str(job_path)})
        replace_bundle(config)
        atomic_json(directory / "status.json", {"ok": True, "backup": config["backup"]})
        return 0
    except Exception as exc:
        atomic_json(directory / "status.json", {"ok": False, "error": str(exc)})
        if parent_exited:
            target = Path(config["target"])
            if target.exists() and bundle_digest(target) == config["originalDigest"]:
                subprocess.run(["/usr/bin/open", str(target)], check=False, capture_output=True)
        return 1


def launch_helper(job: Path) -> None:
    config = json.loads(job.read_text())
    executable = Path(config["staged"]) / "Contents/MacOS/CloningCompanion"
    environment = dict(os.environ, PYINSTALLER_RESET_ENVIRONMENT="1")
    with (job.parent / "helper.log").open("ab") as log:
        process = subprocess.Popen([str(executable), "--apply-local-update", str(job)],
                                   stdout=log, stderr=log, stdin=subprocess.DEVNULL,
                                   start_new_session=True, env=environment)
    deadline = time.monotonic() + 20
    while not (job.parent / "ready.json").exists():
        if process.poll() is not None or time.monotonic() > deadline:
            raise RuntimeError(f"Update helper did not start. App unchanged. See {job.parent / 'helper.log'}")
        time.sleep(.1)
