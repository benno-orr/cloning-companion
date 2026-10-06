# macOS distribution

`scripts/build_macos.sh` builds the self-contained application and DMG. It verifies the nested code signatures before packaging.

## Local artifact

With no signing identity configured, the script applies an ad-hoc signature. This is appropriate for running the app on the build Mac, but Gatekeeper will reject it after public download because it has neither a Developer ID signature nor a notarization ticket.

## Public direct distribution

Requirements:

1. An active Apple Developer Program membership.
2. A `Developer ID Application` certificate installed in the login keychain.
3. A notarytool keychain profile.

Build, sign, submit, and staple:

```bash
SIGN_IDENTITY="Developer ID Application: Your Name (TEAMID)" ./scripts/build_macos.sh
xcrun notarytool store-credentials PlasmidVerifyNotary
NOTARY_PROFILE=PlasmidVerifyNotary ./scripts/notarize_macos.sh
```

The credentials command is interactive; do not put an app-specific password directly in a build script or source file.

## Mac App Store

`app-store.entitlements` is a starting sandbox profile, not a complete App Store release configuration. A store build additionally needs:

- a Mac App Distribution certificate and provisioning profile;
- security-scoped bookmarks for reopening recent files from outside the sandbox;
- App Store Connect metadata, privacy policy, screenshots, and review submission;
- validation of every nested Python framework and extension against App Store signing rules.

The direct-distribution DMG uses `macos.entitlements` and Developer ID signing instead.

