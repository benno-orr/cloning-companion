import base64
import hashlib
import io
import json
import stat
import zipfile
from pathlib import Path

import pytest
from cryptography.exceptions import InvalidSignature
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey

from mac_app import github_updates as updates


@pytest.fixture
def signer(monkeypatch):
    private = Ed25519PrivateKey.generate()
    public = base64.b64encode(private.public_key().public_bytes(serialization.Encoding.Raw, serialization.PublicFormat.Raw)).decode()
    original = updates.verify_manifest
    monkeypatch.setattr(updates, "verify_manifest", lambda envelope: original(envelope, public))
    def sign(**changes):
        manifest = {"schema": 1, "version": "1.4.0", "platform": "macos-arm64",
                    "url": f"https://github.com/{updates.REPOSITORY}/releases/download/v1.4.0/CloningCompanion-1.4.0-arm64.zip",
                    "archive_bytes": 10, "archive_sha256": "a" * 64, "bundle_sha256": "b" * 64}
        manifest.update(changes)
        return {"manifest": manifest, "signature": base64.b64encode(private.sign(updates.canonical(manifest))).decode()}
    return sign


def test_manifest_signature_authenticates_every_field(signer):
    envelope = signer()
    assert updates.verify_manifest(envelope)["version"] == "1.4.0"
    envelope["manifest"]["version"] = "9.0.0"
    with pytest.raises(InvalidSignature):
        updates.verify_manifest(envelope)


@pytest.mark.parametrize("change", [
    {"url": "https://evil.example/app.zip"}, {"platform": "windows"}, {"schema": 2},
    {"archive_bytes": updates.MAX_ARCHIVE + 1}, {"archive_bytes": -1},
    {"bundle_sha256": "x" * 64}, {"version": "1.4.0/../../file"},
])
def test_invalid_signed_metadata_is_still_rejected(signer, change):
    with pytest.raises(ValueError):
        updates.verify_manifest(signer(**change))


def test_remote_check_requires_no_github_credentials(signer, monkeypatch):
    monkeypatch.setattr(updates, "fetch_manifest", signer)
    monkeypatch.setattr(updates.platform, "machine", lambda: "arm64")
    result, envelope = updates.check_update("1.3.0", Path("/Applications/CloningCompanion.app"))
    assert result["available"] and result["canInstall"]
    assert result["channel"] == "GitHub Releases"
    assert not updates.check_update("1.4.0", None)[0]["available"]
    assert not updates.check_update("1.3.0", None)[0]["canInstall"]
    monkeypatch.setattr(updates.platform, "machine", lambda: "x86_64")
    assert not updates.check_update("1.3.0", Path("/Applications/App.app"))[0]["canInstall"]


@pytest.mark.parametrize("url", ["http://github.com/a", "https://evil.example/a", "https://github.com.evil.test/a", "https://user:pass@github.com/a"])
def test_reject_untrusted_downloads_and_redirects(url):
    assert not updates._allowed_url(url)
    with pytest.raises(ValueError):
        updates.open_download(url)
    with pytest.raises(ValueError):
        updates._GitHubRedirect().redirect_request(None, None, 302, "", {}, url)


def make_zip(entries):
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as archive:
        for name, content, mode in entries:
            entry = zipfile.ZipInfo(name)
            entry.external_attr = mode << 16
            archive.writestr(entry, content)
    return output.getvalue()


def test_extract_preserves_executable_and_internal_symlink(tmp_path):
    data = make_zip([
        ("CloningCompanion.app/Contents/MacOS/CloningCompanion", b"executable", stat.S_IFREG | 0o755),
        ("CloningCompanion.app/Contents/link", b"MacOS/CloningCompanion", stat.S_IFLNK | 0o777),
    ])
    path = tmp_path / "app.zip"
    path.write_bytes(data)
    app = updates.extract_app(path, tmp_path / "out")
    assert (app / "Contents/link").read_bytes() == b"executable"
    assert (app / "Contents/MacOS/CloningCompanion").stat().st_mode & 0o111 == 0o111


@pytest.mark.parametrize("entries", [
    [("../../escaped", b"bad", stat.S_IFREG)],
    [("/CloningCompanion.app/escaped", b"bad", stat.S_IFREG)],
    [("Other.app/Contents/bad", b"bad", stat.S_IFREG)],
    [("CloningCompanion.app/Contents/link", b"../../../escape", stat.S_IFLNK)],
    [("CloningCompanion.app/Contents/link", b"/tmp", stat.S_IFLNK)],
    [("CloningCompanion.app/Contents/link", b"elsewhere", stat.S_IFLNK),
     ("CloningCompanion.app/Contents/link/child", b"bad", stat.S_IFREG)],
])
def test_unsafe_archives_are_rejected_before_extraction(tmp_path, entries):
    path = tmp_path / "bad.zip"
    path.write_bytes(make_zip(entries))
    with pytest.raises(ValueError):
        updates.extract_app(path, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_download_verifies_size_hash_then_bundle(signer, tmp_path, monkeypatch):
    data = make_zip([("CloningCompanion.app/Contents/test", b"fixture", stat.S_IFREG)])
    envelope = signer(archive_bytes=len(data), archive_sha256=hashlib.sha256(data).hexdigest())
    monkeypatch.setattr(updates, "open_download", lambda url: io.BytesIO(data))
    checks = []
    monkeypatch.setattr(updates, "verify_bundle", lambda app, manifest: checks.append((app, manifest)))
    manifest = updates.download_update(envelope, tmp_path, "1.3.0")
    assert Path(manifest["bundle"]).is_dir()
    assert checks[0][1]["sha256"] == envelope["manifest"]["bundle_sha256"]
    with pytest.raises(ValueError, match="old release"):
        updates.download_update(envelope, tmp_path, "1.4.0")
    for bad in (data[:-1], data + b"more", b"x" * len(data)):
        monkeypatch.setattr(updates, "open_download", lambda url: io.BytesIO(bad))
        with pytest.raises(ValueError):
            updates.download_update(envelope, tmp_path, "1.3.0")
    assert len(checks) == 1  # No corrupted download reached bundle verification.


def test_manifest_fetch_is_bounded_and_signed(signer, monkeypatch):
    monkeypatch.setattr(updates, "open_download", lambda url: io.BytesIO(json.dumps(signer()).encode()))
    assert updates.fetch_manifest()["manifest"]["version"] == "1.4.0"
    monkeypatch.setattr(updates, "open_download", lambda url: io.BytesIO(b" " * 65537))
    with pytest.raises(ValueError, match="too large"):
        updates.fetch_manifest()
