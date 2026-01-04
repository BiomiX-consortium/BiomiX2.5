#!/usr/bin/env bash
set -e

echo "-----------------------------------------"
echo "Uninstalling BiomiX components..."
echo "-----------------------------------------"

# Get full path to this script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Define paths (relative to installer directory)
R_DIR="$SCRIPT_DIR/R_BiomiX"
PYTHON_DIR="$SCRIPT_DIR/Python_BiomiX"
PYTHON_ENV_DIR="$SCRIPT_DIR/BiomiX-env"
RENV_DIR="$SCRIPT_DIR/renv"
RPROFILE_FILE="$SCRIPT_DIR/.Rprofile"
RBUILD_FILE="$SCRIPT_DIR/_R_build"

# -----------------------------------------
# Helper functions
# -----------------------------------------

delete_if_exists () {
    local TARGET="$1"
    local NAME="$2"

    if [ -d "$TARGET" ]; then
        echo "Deleting $NAME..."
        rm -rf "$TARGET"
        echo "$NAME removed."
    else
        echo "$NAME not found — skipping."
    fi
}

delete_file_if_exists () {
    local FILE="$1"
    local NAME="$2"

    if [ -f "$FILE" ]; then
        echo "Deleting $NAME..."
        rm -f "$FILE"
        echo "$NAME removed."
    else
        echo "$NAME not found — skipping."
    fi
}

# -----------------------------------------
# Uninstallation steps
# -----------------------------------------

delete_if_exists "$R_DIR" "R_BiomiX"
delete_if_exists "$PYTHON_ENV_DIR" "Python BiomiX-env"
delete_if_exists "$PYTHON_DIR" "Python_BiomiX"
delete_if_exists "$RENV_DIR" "renv environment"
delete_file_if_exists "$RPROFILE_FILE" ".Rprofile file"
delete_if_exists "$RBUILD_FILE" "_R_build folder"

echo
echo "Uninstallation complete."
echo "Press ENTER to exit."
read
exit 0

