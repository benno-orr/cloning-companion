"""Package and sign a release locally. This script never uploads anything."""
import argparse
import base64
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from mac_app.github_updates import REPOSITORY, canonical, verify_manifest
from mac_app.local_updates import bundle_info, bundle_digest


def signing_key():
    directory = Path.home() / "Library/Application Support/CloningCompanion/Signing"
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    path = directory / "release-ed25519.pem"
    if not path.exists():
        key = Ed25519PrivateKey.generate()
        payload = key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8, serialization.NoEncryption())
        with os.fdopen(os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600), "wb") as handle:
            handle.write(payload)
    return serialization.load_pem_private_key(path.read_bytes(), password=None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("app", nargs="?", type=Path)
    parser.add_argument("--public-key", action="store_true")
    args = parser.parse_args()
    key = signing_key()
    public = base64.b64encode(key.public_key().public_bytes(serialization.Encoding.Raw, serialization.PublicFormat.Raw)).decode()
    if args.public_key:
        print(public)
        return
    if args.app is None:
        parser.error("app path required")
    version = bundle_info(args.app)["CFBundleShortVersionString"]
    directory = args.app.parent / f"github-release-{version}"
    directory.mkdir(exist_ok=True)
    archive = directory / f"CloningCompanion-{version}-arm64.zip"
    if archive.exists():
        raise ValueError("Release archive already exists; use a new version or move the existing archive first")
    subprocess.run(["/usr/bin/codesign", "--verify", "--deep", "--strict", str(args.app)], check=True)
    subprocess.run(["/usr/bin/ditto", "-c", "-k", "--norsrc", "--noextattr", "--keepParent", str(args.app), str(archive)], check=True)
    manifest = {"schema": 1, "version": version, "platform": "macos-arm64", "archive_bytes": archive.stat().st_size,
                "url": f"https://github.com/{REPOSITORY}/releases/download/v{version}/{archive.name}",
                "archive_sha256": hashlib.file_digest(archive.open("rb"), "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(archive.read_bytes()).hexdigest(),
                "bundle_sha256": bundle_digest(args.app)}
    envelope = {"manifest": manifest, "signature": base64.b64encode(key.sign(canonical(manifest))).decode()}
    verify_manifest(envelope)
    (directory / "latest.json").write_text(json.dumps(envelope, indent=2) + "\n")
    print(directory)


if __name__ == "__main__":
    main()
