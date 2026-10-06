#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="${0:A:h}"
PROJECT_DIR="${SCRIPT_DIR:h}"
PYTHON_BIN="${PYTHON_BIN:-${PROJECT_DIR}/.venv/bin/python}"
APP_PATH="${PROJECT_DIR}/dist/CloningCompanion.app"
DMG_STAGE="${PROJECT_DIR}/build/dmg-root-1.5.1"
DMG_PATH="${PROJECT_DIR}/dist/CloningCompanion-1.5.1.dmg"
SIGN_IDENTITY="${SIGN_IDENTITY:-}"

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "Python environment not found at ${PYTHON_BIN}"
  echo "Create it and install requirements-macos.txt first."
  exit 1
fi

cd "${PROJECT_DIR}"
"${PYTHON_BIN}" -m PyInstaller --noconfirm --clean PlasmidVerify.spec

xattr -cr "${APP_PATH}"
if [[ -n "${SIGN_IDENTITY}" ]]; then
  codesign --force --deep --options runtime --timestamp \
    --entitlements "${PROJECT_DIR}/packaging/macos.entitlements" \
    --sign "${SIGN_IDENTITY}" "${APP_PATH}"
else
  echo "SIGN_IDENTITY is unset; applying an ad-hoc signature for local testing."
  codesign --force --deep --sign - "${APP_PATH}"
fi

codesign --verify --deep --strict --verbose=2 "${APP_PATH}"
"${PYTHON_BIN}" "${PROJECT_DIR}/scripts/publish_local_update.py" "${APP_PATH}"

# Day-to-day development updates no longer need a disk image.
if [[ "${BUILD_DMG:-1}" == "0" ]]; then
  echo "Built and published ${APP_PATH} for in-app updating"
  exit 0
fi

rm -rf "${DMG_STAGE}"
mkdir -p "${DMG_STAGE}"
ditto "${APP_PATH}" "${DMG_STAGE}/CloningCompanion.app"
ln -s /Applications "${DMG_STAGE}/Applications"
rm -f "${DMG_PATH}"
hdiutil create -volname "CloningCompanion" -srcfolder "${DMG_STAGE}" \
  -ov -format UDZO "${DMG_PATH}"

if [[ -n "${SIGN_IDENTITY}" ]]; then
  codesign --force --timestamp --sign "${SIGN_IDENTITY}" "${DMG_PATH}"
fi

echo "Built ${APP_PATH}"
echo "Built ${DMG_PATH}"
