#!/usr/bin/env bash
set -e

########################################
# Settings
########################################

R_VERSION="4.4.1"
CRAN_MIRROR="https://cloud.r-project.org"

# Directory dello script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

########################################
# Detect OS
########################################

if ! command -v lsb_release >/dev/null 2>&1; then
    echo "lsb_release not found. Unsupported system."
    read -p "Press Enter to exit..."
    exit 1
fi

DISTRO="$(lsb_release -is)"
CODENAME="$(lsb_release -cs)"

if [[ "$DISTRO" != "Ubuntu" && "$DISTRO" != "Debian" ]]; then
    echo "Unsupported distribution: $DISTRO"
    read -p "Press Enter to exit..."
    exit 1
fi

########################################
# Clean previous R packages (if any)
########################################

echo "Cleaning any previous R installation..."
sudo apt purge -y r-base r-recommended '^r-cran-.*' || true
sudo apt autoremove -y || true

########################################
# Install R core only
########################################

echo "Installing R $R_VERSION (r-base-core only)..."

sudo apt update

sudo apt install -y --no-install-recommends \
  r-base-core="$R_VERSION-3.2204.0"

# Lock the version to prevent upgrades
sudo apt-mark hold r-base-core

sudo apt install -y libblas-dev liblapack-dev

sudo apt install -y zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

sudo apt install -y libudunits2-dev

sudo apt install -y \
  libgdal-dev \
  libgeos-dev \
  libproj-dev \
  gdal-bin \
  proj-bin \
  proj-data


########################################
# Verify installation
########################################

INSTALLED_VERSION="$(R --version | head -n1 | awk '{print $3}')"
if [[ "$INSTALLED_VERSION" != "$R_VERSION" ]]; then
    echo "ERROR: Installed R version $INSTALLED_VERSION (expected $R_VERSION)"
    read -p "Press Enter to exit..."
    exit 1
fi

echo "R $INSTALLED_VERSION successfully installed!"
echo "R binary: $(command -v R)"
echo "Rscript: $(command -v Rscript)"

########################################
# Install renv
########################################

echo "Installing renv package..."
Rscript -e "install.packages('renv', repos='$CRAN_MIRROR')"

########################################
# Initialize renv in the project folder
########################################

echo "Initializing renv in project folder..."
Rscript -e "renv::init(bare = TRUE)"

########################################
# Done
########################################

echo
echo "R environment setup complete!"
read -p "Press Enter to continue..."

