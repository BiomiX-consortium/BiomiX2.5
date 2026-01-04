#!/usr/bin/env bash
set -e

Rscript --version


########################################
# Deactivate flag for Rdisop package installation
########################################

BUILD_DIR="$(pwd)/_R_build"
mkdir -p "$BUILD_DIR"

cat <<EOF > "$BUILD_DIR/Makevars"
CFLAGS = -g -O2 -Wno-error=format-security
CXXFLAGS = -g -O2 -Wno-error=format-security
CXX17FLAGS = -g -O2 -Wno-error=format-security
EOF

export R_MAKEVARS_USER="$BUILD_DIR/Makevars"



########################################
# Get full path to this script's folder
########################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

########################################
# Detect Rscript
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
"$R_EXE" "$SCRIPT_DIR/MODULE_LINUX_R_packages.r"

########################################
# Pause
########################################

read -p "Press Enter to continue..."

