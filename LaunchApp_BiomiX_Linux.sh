#!/bin/bash
# =============================================================================
# launch_biomix.sh
# BiomiX launcher script for Linux
#
# Launches the Shiny app using the renv environment.
# =============================================================================

# Get the directory where this script lives (= BiomiX root directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Open the browser automatically after a short delay (waits for Shiny to start)
sleep 2 && xdg-open http://127.0.0.1:3838 &

# Launch the Shiny app using the renv R installation
"${SCRIPT_DIR}/_INSTALL/R_BiomiX/bin/Rscript" -e \
  "renv::load('${SCRIPT_DIR}/_INSTALL'); setwd('${SCRIPT_DIR}'); shiny::runApp('.', port=3838, launch.browser=FALSE)"