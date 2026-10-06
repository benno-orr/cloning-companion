#!/bin/zsh
set -euo pipefail

SCRIPT_DIR="${0:A:h}"
PROJECT_DIR="${SCRIPT_DIR:h}"
DMG_PATH="${1:-${PROJECT_DIR}/dist/CloningCompanion-1.4.0.dmg}"
NOTARY_PROFILE="${NOTARY_PROFILE:-}"

if [[ -z "${NOTARY_PROFILE}" ]]; then
  echo "Set NOTARY_PROFILE to a notarytool keychain profile name."
  echo "Create one with: xcrun notarytool store-credentials PROFILE_NAME"
  exit 1
fi

if [[ ! -f "${DMG_PATH}" ]]; then
  echo "DMG not found: ${DMG_PATH}"
  exit 1
fi

xcrun notarytool submit "${DMG_PATH}" --keychain-profile "${NOTARY_PROFILE}" --wait
xcrun stapler staple "${DMG_PATH}"
xcrun stapler validate "${DMG_PATH}"
spctl --assess --type open --context context:primary-signature --verbose=2 "${DMG_PATH}"

echo "Notarized and stapled ${DMG_PATH}"
