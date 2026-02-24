#!/bin/bash
set -e

# System-wide conda prefix. Override with: CONDA_BASE=/other bash BiomiX_INSTALL_LINUX.sh
CONDA_BASE="${CONDA_BASE:-/opt/conda}"

source "$CONDA_BASE/etc/profile.d/conda.sh"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

conda env list | grep -q '^biomix ' && conda env remove -n biomix --yes

# Create the full environment in one solver pass
mamba env create -f "$SCRIPT_DIR/biomix_env.yml"

conda activate biomix

# Install the small set of custom R packages not available on conda
Rscript "$SCRIPT_DIR/MODULE_LINUX_custom_r.r"
