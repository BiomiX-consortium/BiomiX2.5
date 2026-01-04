#!/usr/bin/env bash
set -e


########################################
# Deactivate flag for sp package installation
########################################
INSTALL_DIR=$(pwd)/_R_build
mkdir -p "$INSTALL_DIR"

cat <<EOF > "$INSTALL_DIR/Makevars"
CFLAGS = -g -O2 -Wno-error=format-security
CXXFLAGS = -g -O2 -Wno-error=format-security
EOF

export R_MAKEVARS_USER="$INSTALL_DIR/Makevars"

########################################
# Get full path to this script's folder
########################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

########################################
# Detect Rscript (Linux equivalent)
########################################

R_EXE=""

if command -v Rscript >/dev/null 2>&1; then
    R_EXE="Rscript"
else
    echo "Rscript not found on system"
    read -p "Press Enter to exit..."
    exit 1
fi

########################################
# Test R is working
########################################

"$R_EXE" -e "cat('R_BiomiX Rscript detected')" || echo "Failed to run R_BiomiX Rscript"

########################################
# Run the script
########################################

echo "Running install_from_github.R..."
"$R_EXE" "$SCRIPT_DIR/MODULE_LINUX_R_Extension_SNF_NEMO.r"

########################################
# Pause
########################################

read -p "Press Enter to continue..."

