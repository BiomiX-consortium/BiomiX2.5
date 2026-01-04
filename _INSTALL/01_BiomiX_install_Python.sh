#!/usr/bin/env bash
set -e

########################################
# Settings
########################################

PYTHON_VERSION="3.9.13"

# Directory dello script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Install dir (manteniamo la stessa struttura logica)
INSTALL_DIR="$SCRIPT_DIR/Python_BiomiX"

# Environment name and path
ENV_NAME="BiomiX-env"
ENV_DIR="$SCRIPT_DIR/$ENV_NAME"

########################################
# Detect system architecture
########################################

ARCHITECTURE="x86"
UNAME_ARCH="$(uname -m)"

if [[ "$UNAME_ARCH" == "x86_64" || "$UNAME_ARCH" == "amd64" ]]; then
    ARCHITECTURE="amd64"
fi

echo "$ARCHITECTURE"

########################################
# Detect Python executable
########################################

PYTHON_BIN=""

if command -v python3.9 >/dev/null 2>&1; then
    PYTHON_BIN="python3.9"
elif command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
else
    echo "Python 3 not found on system."
    read -p "Press Enter to exit..."
    exit 1
fi

########################################
# Verify Python version
########################################

PY_VERSION_STR="$($PYTHON_BIN --version 2>&1 | awk '{print $2}')"

if [[ "$PY_VERSION_STR" != 3.9.* ]]; then
    echo "WARNING: Detected Python version $PY_VERSION_STR"
    echo "Expected Python $PYTHON_VERSION"
fi

########################################
# Create virtual environment
########################################

echo "Creating virtual environment: $ENV_NAME"
"$PYTHON_BIN" -m venv "$ENV_DIR"

########################################
# Activate the environment
########################################

# shellcheck source=/dev/null
source "$ENV_DIR/bin/activate"

########################################
# Upgrade pip and install packages
########################################

python -m pip install --upgrade pip setuptools wheel

########################################
# Install packages
########################################

python -m pip install mamba

python -m pip install \
    PyQt5==5.12.3 \
    pandas==2.2.3 \
    scikit-learn==1.6.1 \
    xlrd==2.0.1 \
    openpyxl==3.1.5 \
    mofapy2==0.7.1 \
    numpy==1.26.4

########################################
# Completion message
########################################

echo
echo "Python Environment setup complete!"
echo "To activate it in the future, run:"
echo "source \"$ENV_DIR/bin/activate\""
read -p "Press Enter to continue..."

