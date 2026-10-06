# -*- mode: python ; coding: utf-8 -*-
from pathlib import Path

from PyInstaller.utils.hooks import collect_all


root = Path(SPECPATH)
webview_datas, webview_binaries, webview_hidden = collect_all("webview")

a = Analysis(
    [str(root / "mac_app" / "main.py")],
    pathex=[str(root)],
    binaries=webview_binaries,
    datas=webview_datas + [
        (str(root / "mac_app" / "index.html"), "mac_app"),
        (str(root / "assets" / "CloningCompanion-logo.png"), "mac_app/assets"),
    ],
    hiddenimports=webview_hidden + ["webview.platforms.cocoa"],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["tkinter", "pytest", "streamlit", "pandas"],
    noarchive=False,
    optimize=1,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    name="CloningCompanion",
    debug=False,
    exclude_binaries=True,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch="arm64",
    codesign_identity=None,
    entitlements_file=str(root / "packaging" / "macos.entitlements"),
)

collection = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="CloningCompanion",
)

app = BUNDLE(
    collection,
    name="CloningCompanion.app",
    icon=str(root / "assets" / "CloningCompanion.icns"),
    bundle_identifier="com.bennoorr.cloningcompanion",
    version="1.5.1",
    info_plist={
        "CFBundleDisplayName": "CloningCompanion",
        "CFBundleName": "CloningCompanion",
        "CFBundleShortVersionString": "1.5.1",
        "CFBundleVersion": "33",
        "NSHumanReadableCopyright": "Copyright © 2026 Benno Orr. All rights reserved.",
        "NSHighResolutionCapable": True,
        "LSMinimumSystemVersion": "12.0",
        "NSPrincipalClass": "NSApplication",
        "CFBundleDocumentTypes": [
            {
                "CFBundleTypeName": "CloningCompanion Project",
                "CFBundleTypeExtensions": ["plasmidverify"],
                "CFBundleTypeRole": "Editor",
                "LSHandlerRank": "Owner",
            }
        ],
    },
)
