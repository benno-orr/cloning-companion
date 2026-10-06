import json
import plistlib
import shutil
import subprocess
from pathlib import Path

import pytest

from mac_app import local_updates as updates


def app(path, version="1.2.4"):
    (path / "Contents/MacOS").mkdir(parents=True)
    (path / "Contents/Info.plist").write_bytes(plistlib.dumps({
        "CFBundleIdentifier": updates.BUNDLE_ID, "CFBundleExecutable": "CloningCompanion",
        "CFBundleShortVersionString": version}))
    (path / "Contents/MacOS/CloningCompanion").write_text("fixture " + version)
    return path


@pytest.fixture
def environment(tmp_path, monkeypatch):
    directory = tmp_path / "updates"
    directory.mkdir()
    target = app(tmp_path / "CloningCompanion.app")
    source = app(tmp_path / "build" / "CloningCompanion.app", "1.3.0")
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        if command[0] == "/usr/bin/ditto":
            shutil.copytree(command[1], command[2], symlinks=True)
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(updates.subprocess, "run", run)
    monkeypatch.setattr(updates, "update_directory", lambda: directory)
    updates.publish(source, directory)
    return directory, target, source, calls


def test_update_discovery_versions_and_preview_restriction(environment):
    directory, target, _, _ = environment
    assert updates.check_update("1.2.4", directory, target)["canInstall"]
    assert not updates.check_update("1.2.4", directory, None)["canInstall"]
    assert not updates.check_update("1.3.0", directory, target)["available"]
    assert not updates.check_update("1.10.0", directory, target)["available"]
    assert not updates.check_update("1.2.4", directory, Path('/Volumes/Test/App.app'))["canInstall"]
    with pytest.raises(ValueError):
        updates.version_key("1.3.0; shell command")


def test_successful_replacement_preserves_backup_and_unsaved_inputs(environment):
    directory, target, _, _ = environment
    original = updates.bundle_digest(target)
    session = {"state": {"inserts": [{"sequence": "ATG"}]}, "path": None}
    job = updates.stage_update(target, directory, "1.2.4", session)
    assert updates.bundle_digest(target) == original  # Staging never edits the installed app.
    assert json.loads((job.parent / "recovery.plasmidverify").read_text())["state"] == session["state"]
    config = json.loads(job.read_text())
    launched = []
    updates.replace_bundle(config, launch=launched.append)
    assert updates.bundle_info(target)["CFBundleShortVersionString"] == "1.3.0"
    assert updates.bundle_digest(Path(config["backup"])) == original
    assert launched == [target]


def test_launch_failure_rolls_back_and_reopens_old_app(environment):
    directory, target, _, _ = environment
    original = updates.bundle_digest(target)
    config = json.loads(updates.stage_update(target, directory, "1.2.4", {"state": {}}).read_text())
    attempts = []
    def launch(path):
        attempts.append(updates.bundle_info(path)["CFBundleShortVersionString"])
        if len(attempts) == 1:
            raise RuntimeError("launch failed")
    with pytest.raises(RuntimeError, match="launch failed"):
        updates.replace_bundle(config, launch)
    assert updates.bundle_digest(target) == original
    assert attempts == ["1.3.0", "1.2.4"]


def test_corrupt_publication_does_not_touch_installed_app(environment):
    directory, target, source, _ = environment
    original = updates.bundle_digest(target)
    (source / "Contents/MacOS/CloningCompanion").write_text("modified")
    with pytest.raises(ValueError, match="checksum"):
        updates.stage_update(target, directory, "1.2.4", {"state": {}})
    assert updates.bundle_digest(target) == original


def test_recheck_staged_and_installed_bundles_before_replacement(environment):
    directory, target, _, _ = environment
    config = json.loads(updates.stage_update(target, directory, "1.2.4", {"state": {}}).read_text())
    (Path(config["staged"]) / "Contents/MacOS/CloningCompanion").write_text("modified")
    with pytest.raises(ValueError, match="checksum"):
        updates.replace_bundle(config)
    assert updates.bundle_info(target)["CFBundleShortVersionString"] == "1.2.4"
    (target / "Contents/MacOS/CloningCompanion").write_text("changed externally")
    with pytest.raises(ValueError, match="changed while"):
        updates.replace_bundle(config)


def test_reject_other_app_and_external_symlinks(environment):
    _, target, _, _ = environment
    (target / "Contents/outside").symlink_to(target.parent)
    with pytest.raises(ValueError):
        updates.bundle_digest(target)
    info_path = target / "Contents/Info.plist"
    info = plistlib.loads(info_path.read_bytes())
    info["CFBundleIdentifier"] = "someone.else"
    info_path.write_bytes(plistlib.dumps(info))
    with pytest.raises(ValueError, match="not a CloningCompanion"):
        updates.bundle_info(target)


def test_failed_signature_stops_staging(environment, monkeypatch):
    directory, target, _, _ = environment
    original = updates.bundle_digest(target)
    real_run = updates.subprocess.run
    def run(command, **kwargs):
        if command[0] == "/usr/bin/codesign":
            raise subprocess.CalledProcessError(1, command)
        return real_run(command, **kwargs)
    monkeypatch.setattr(updates.subprocess, "run", run)
    with pytest.raises(subprocess.CalledProcessError):
        updates.stage_update(target, directory, "1.2.4", {"state": {}})
    assert updates.bundle_digest(target) == original


def test_helper_waits_for_exit_then_replaces(environment, monkeypatch):
    directory, target, _, _ = environment
    job = updates.stage_update(target, directory, "1.2.4", {"state": {"inserts": []}})
    events = []
    def kill(pid, signal):
        events.append("poll")
        if len(events) > 1:
            raise ProcessLookupError()
    monkeypatch.setattr(updates.os, "kill", kill)
    monkeypatch.setattr(updates.time, "sleep", lambda _: None)
    assert updates.apply_update_job(job) == 0
    assert len(events) == 2
    assert json.loads((job.parent / "status.json").read_text())["ok"]
    assert json.loads((directory / "pending-session.json").read_text())["job"] == str(job)


def test_session_recovery_is_retained_until_ui_acknowledges(environment, monkeypatch):
    from mac_app.main import NativeAPI
    directory, target, _, _ = environment
    job = updates.stage_update(target, directory, "1.2.4", {"state": {"inserts": []}, "path": None})
    updates.atomic_json(directory / "pending-session.json", {"job": str(job)})
    api = NativeAPI()
    monkeypatch.setattr(api, "recent_projects", lambda: [])
    assert api.app_info()["updateRecovery"]["state"] == {"inserts": []}
    assert (directory / "pending-session.json").exists()
    api.acknowledge_update_recovery()
    assert "updateRecovery" not in api.app_info()
    assert (job.parent / "recovery.plasmidverify").exists()


@pytest.mark.parametrize("fail_launch", [False, True])
def test_contents_swap_for_writable_app_in_protected_parent(environment, monkeypatch, fail_launch):
    directory, target, _, _ = environment
    original = updates.bundle_digest(target)
    real_mkdtemp = updates.tempfile.mkdtemp
    def mkdtemp(*args, **kwargs):
        if kwargs.get("dir") == target.parent:
            raise PermissionError("Protected Applications folder")
        return real_mkdtemp(*args, **kwargs)
    monkeypatch.setattr(updates.tempfile, "mkdtemp", mkdtemp)
    job = updates.stage_update(target, directory, "1.2.4", {"state": {}})
    config = json.loads(job.read_text())
    assert config["mode"] == "contents"
    assert Path(config["backup"]).is_relative_to(directory)
    attempts = []
    def launch(path):
        attempts.append(path)
        if fail_launch and len(attempts) == 1:
            raise RuntimeError("failed")
    if fail_launch:
        with pytest.raises(RuntimeError):
            updates.replace_bundle(config, launch)
        assert updates.bundle_digest(target) == original
    else:
        updates.replace_bundle(config, launch)
        assert updates.bundle_info(target)["CFBundleShortVersionString"] == "1.3.0"
        assert updates.bundle_digest(Path(config["backup"])) == original
    assert {p.name for p in target.iterdir()} == {"Contents"}
