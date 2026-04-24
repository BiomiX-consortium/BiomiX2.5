#!/bin/bash
# build_pkg.sh — Developer script to build BiomiX2.5.pkg
# Run from the repository root: bash _INSTALL/mac/build_pkg.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
MAC_DIR="$SCRIPT_DIR"
PAYLOAD_DIR="$MAC_DIR/payload"
APP_CONTENTS="$PAYLOAD_DIR/Applications/BiomiX.app/Contents"
RESOURCES_DIR="$APP_CONTENTS/Resources"
BIOMIX_SRC="$RESOURCES_DIR/BiomiX"
OUTPUT_PKG="$REPO_ROOT/BiomiX2.5.pkg"

echo "=== BiomiX2.5 macOS Package Builder ==="
echo "Repo root: $REPO_ROOT"

# --- Clean previous payload resources (may be root-owned if a previous install ran postinstall) ---
echo "Cleaning previous build artifacts..."
sudo rm -rf "$BIOMIX_SRC"
sudo rm -f "$RESOURCES_DIR/biomix_env_mac.yml" \
           "$RESOURCES_DIR/biomix_env.yml" \
           "$RESOURCES_DIR/MODULE_LINUX_custom_r.r" \
           "$RESOURCES_DIR/MODULE_LINUX_R_Extension_SNF_NEMO.r" \
           "$RESOURCES_DIR/MODULE_LINUX_R_Extension_MintTea.r" \
           "$RESOURCES_DIR/MODULE_MAC_custom_r.r" \
           "$RESOURCES_DIR/BiomiX.icns"
# Restore ownership of the payload tree so subsequent commands run without sudo
sudo chown -R "$(whoami)" "$PAYLOAD_DIR"

mkdir -p "$BIOMIX_SRC"

# --- Copy BiomiX source files (exclude Windows-only, mac build dir, git, etc.) ---
echo "Copying BiomiX source files..."
rsync -a \
    --exclude '_INSTALL/' \
    --exclude '__pycache__/' \
    --exclude '.git/' \
    --exclude '.DS_Store' \
    --exclude 'Rplots.pdf' \
    --exclude 'BiomiX2.5.pkg' \
    --exclude 'LaunchApp_BiomiX_Windows.bat' \
    --exclude 'MODULE_WINDOWS.py' \
    --exclude 'directory.txt' \
    --exclude 'directory_out.txt' \
    --exclude 'COMBINED_COMMANDS.json' \
    "$REPO_ROOT/" "$BIOMIX_SRC/"

# --- Copy install resources ---
echo "Copying install resources..."
cp "$MAC_DIR/BiomiX.icns"                                      "$RESOURCES_DIR/"
cp "$REPO_ROOT/_INSTALL/biomix_env_mac.yml"                    "$RESOURCES_DIR/"
cp "$REPO_ROOT/_INSTALL/MODULE_LINUX_custom_r.r"               "$RESOURCES_DIR/"
cp "$REPO_ROOT/_INSTALL/MODULE_LINUX_R_Extension_SNF_NEMO.r"   "$RESOURCES_DIR/"
cp "$REPO_ROOT/_INSTALL/MODULE_LINUX_R_Extension_MintTea.r"    "$RESOURCES_DIR/"
cp "$REPO_ROOT/_INSTALL/MODULE_MAC_custom_r.r"                 "$RESOURCES_DIR/"

# --- Generate component plist with relocation disabled ---
COMPONENT_PLIST="$MAC_DIR/component.plist"
cat > "$COMPONENT_PLIST" <<'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<array>
    <dict>
        <key>BundleHasStrictIdentifier</key>
        <true/>
        <key>BundleIsRelocatable</key>
        <false/>
        <key>BundleIsVersionChecked</key>
        <true/>
        <key>BundleOverwriteAction</key>
        <string>upgrade</string>
        <key>RootRelativeBundlePath</key>
        <string>Applications/BiomiX.app</string>
    </dict>
</array>
</plist>
EOF

# --- Build the .pkg ---
echo "Building $OUTPUT_PKG..."
pkgbuild \
    --root "$PAYLOAD_DIR" \
    --scripts "$MAC_DIR/scripts" \
    --component-plist "$COMPONENT_PLIST" \
    --identifier com.biomix.BiomiX \
    --version 2.5 \
    --install-location / \
    "$OUTPUT_PKG"

echo ""
echo "=== Done ==="
echo "Package: $OUTPUT_PKG"
echo ""
echo "To test locally (no signing required for lab distribution):"
echo "  sudo installer -pkg BiomiX2.5.pkg -target /"
