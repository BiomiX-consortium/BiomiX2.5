#!/bin/bash
CONDA_BASE="${CONDA_BASE:-/opt/conda}"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate biomix
python3 MODULE_LINUX.py
