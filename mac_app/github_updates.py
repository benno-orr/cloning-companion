"""Public GitHub Releases updates, authenticated by a pinned Ed25519 key."""
from __future__ import annotations

import base64
import hashlib
import json
import os
import platform
import posixpath
import ssl
import stat
import tempfile
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

import certifi
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey

from mac_app.local_updates import version_key, verify_bundle
from mac_app.update_identity import PUBLIC_KEY

REPOSITORY = "benno-orr/cloning-companion"
FEED_URL = f"https://github.com/{REPOSITORY}/releases/latest/download/latest.json"
MAX_ARCHIVE = 1024 * 1024 * 1024
MAX_EXPANDED = 2 * MAX_ARCHIVE


def canonical(data: dict) -> bytes:
    return json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def verify_manifest(envelope: dict, public_key: str = PUBLIC_KEY) -> dict:
    manifest = envelope["manifest"]
    Ed25519PublicKey.from_public_bytes(base64.b64decode(public_key, validate=True)).verify(
        base64.b64decode(envelope["signature"], validate=True), canonical(manifest))
    if manifest.get("schema") != 1 or manifest.get("platform") != "macos-arm64":
        raise ValueError("Unsupported update format or platform")
    version = manifest["version"]
    version_key(version)
    expected = f"https://github.com/{REPOSITORY}/releases/download/v{version}/CloningCompanion-{version}-arm64.zip"
    if manifest["url"] != expected:
        raise ValueError("Update asset is not from the configured repository/version")
    if type(manifest["archive_bytes"]) is not int or not 0 < manifest["archive_bytes"] <= MAX_ARCHIVE:
        raise ValueError("Invalid update download size")
    for key in ("archive_sha256", "bundle_sha256"):
        if len(manifest[key]) != 64 or any(c not in "0123456789abcdef" for c in manifest[key]):
            raise ValueError("Invalid update checksum")
    return manifest


def _allowed_url(url: str) -> bool:
    parts = urlsplit(url)
    return (parts.scheme == "https" and not parts.username and not parts.password and
            parts.port in (None, 443) and parts.hostname in {
                "github.com", "release-assets.githubusercontent.com", "objects.githubusercontent.com"})


class _GitHubRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if not _allowed_url(newurl):
            raise ValueError("Refusing an update redirect outside GitHub HTTPS")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def open_download(url: str):
    if not _allowed_url(url):
        raise ValueError("Updates require GitHub HTTPS")
    opener = urllib.request.build_opener(_GitHubRedirect(), urllib.request.HTTPSHandler(
        context=ssl.create_default_context(cafile=certifi.where())))
    return opener.open(urllib.request.Request(url, headers={"User-Agent": "CloningCompanion-Updater", "Cache-Control": "no-cache"}), timeout=30)


def fetch_manifest() -> dict:
    with open_download(FEED_URL) as response:
        data = response.read(65537)
    if len(data) > 65536:
        raise ValueError("Update metadata is too large")
    envelope = json.loads(data)
    verify_manifest(envelope)
    return envelope


def check_update(current: str, target: Path | None) -> tuple[dict, dict]:
    envelope = fetch_manifest()
    manifest = verify_manifest(envelope)
    available = version_key(manifest["version"]) > version_key(current)
    reason = None
    if platform.machine() not in ("arm64", "aarch64"):
        reason = "This release supports Apple silicon Macs only."
    elif target is None:
        reason = "Open the desktop app to install updates."
    elif str(target).startswith("/Volumes/") or "AppTranslocation" in target.parts:
        reason = "Move the app to Applications before updating."
    return {"ok": True, "available": available, "version": manifest["version"], "canInstall": reason is None,
            "channel": "GitHub Releases", "message": reason if available and reason else
            f"Version {manifest['version']} is ready from GitHub." if available else
            f"You’re up to date with GitHub releases ({current})."}, envelope


def extract_app(archive: Path, destination: Path) -> Path:
    """Extract only one .app; reject traversal, symlink escapes and zip bombs."""
    with zipfile.ZipFile(archive) as package:
        members = package.infolist()
        if len(members) > 100000 or sum(item.file_size for item in members) > MAX_EXPANDED:
            raise ValueError("Update expands beyond the allowed limit")
        seen, links = set(), {}
        for item in members:
            name = item.filename.rstrip("/")
            parts = PurePosixPath(name).parts
            if (not parts or parts[0] != "CloningCompanion.app" or ".." in parts or
                    "\\" in name or name.startswith("/") or name in seen):
                raise ValueError("Unsafe or duplicate archive path")
            seen.add(name)
            mode = item.external_attr >> 16
            if stat.S_ISLNK(mode):
                if item.file_size > 4096:
                    raise ValueError("Invalid symlink")
                link = package.read(item).decode("utf-8")
                resolved = posixpath.normpath(posixpath.join(posixpath.dirname(name), link))
                if link.startswith("/") or not resolved.startswith("CloningCompanion.app/"):
                    raise ValueError("Update symlink escapes the app")
                links[name] = link
            elif stat.S_IFMT(mode) not in (0, stat.S_IFREG, stat.S_IFDIR):
                raise ValueError("Unsupported archive entry type")
        for name in seen:
            if any(str(parent) in links for parent in PurePosixPath(name).parents):
                raise ValueError("Archive entry traverses a symlink")
        for item in members:
            name = item.filename.rstrip("/")
            path = destination / name
            if name in links:
                continue
            if item.is_dir():
                path.mkdir(parents=True, exist_ok=True)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                with package.open(item) as source, path.open("xb") as output:
                    for block in iter(lambda: source.read(1024 * 1024), b""):
                        output.write(block)
                path.chmod(0o644 | ((item.external_attr >> 16) & 0o111))
        for name, link in links.items():
            path = destination / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.symlink_to(link)
    return destination / "CloningCompanion.app"


def download_update(envelope: dict, directory: Path, current: str) -> dict:
    manifest = verify_manifest(envelope)
    if version_key(manifest["version"]) <= version_key(current):
        raise ValueError("Refusing to install an old release")
    folder = Path(tempfile.mkdtemp(prefix="github-", dir=directory))
    archive = folder / "update.zip"
    digest, received = hashlib.sha256(), 0
    with open_download(manifest["url"]) as response, archive.open("xb") as output:
        for block in iter(lambda: response.read(1024 * 1024), b""):
            received += len(block)
            if received > manifest["archive_bytes"]:
                raise ValueError("Update download exceeds its signed size")
            digest.update(block)
            output.write(block)
    if received != manifest["archive_bytes"] or digest.hexdigest() != manifest["archive_sha256"]:
        raise ValueError("Update download is incomplete or corrupted")
    app = extract_app(archive, folder / "expanded")
    local_manifest = {"version": manifest["version"], "bundle": str(app), "sha256": manifest["bundle_sha256"], "channel": "github"}
    verify_bundle(app, local_manifest)
    return local_manifest
