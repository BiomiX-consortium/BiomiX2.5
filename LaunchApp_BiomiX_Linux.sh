#!/bin/bash
CONDA_BASE="${CONDA_BASE:-/opt/conda}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate Biomix2.5
python3 MODULE_LINUX.py
