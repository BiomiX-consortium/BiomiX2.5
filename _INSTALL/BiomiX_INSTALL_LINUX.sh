#!/bin/bash
set -e

# System-wide conda prefix. Override with: CONDA_BASE=/other bash BiomiX_INSTALL_LINUX.sh
CONDA_BASE="${CONDA_BASE:-/opt/conda}"

source "$CONDA_BASE/etc/profile.d/conda.sh"

conda env list | grep -q '^biomix ' && conda env remove -n biomix --yes

conda create -n biomix python=3.9 -y
conda activate biomix
conda config --add channels conda-forge
conda install mamba=2.0.5 -y
mamba install -c conda-forge r-base=4.4.1 -y
mamba install -c conda-forge r-systemfonts=1.1.0 -y
mamba install -c conda-forge r-rcpp=1.0.13 -y
python3 MODULE_LINUX.py
