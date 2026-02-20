# BiomiX – MintTea R extension installer
# Installs the MintTea package from GitHub.
# Reference: Muller et al., Nature Communications 2024
# https://github.com/efratmuller/MintTea
#
# NOTE: This step requires an active internet connection.

# # MANUAL INPUT (for testing)
# args = list("CASE", "CTRL", "/path/to/BiomiX2.5")
# directory <- args[[3]]
# setwd(paste(directory, "_INSTALL", sep = "/"))
# renv::load(paste(directory, "_INSTALL", sep = "/"))


# ---- Activate renv environment ----
cat("Activating renv environment...\n")
renv::activate()

cat("Installing into library:\n")
print(.libPaths())

chooseCRANmirror(48, ind = TRUE)
options(Ncpus = 6)


# ---- Install MintTea dependencies not yet bundled ----
# MintTea depends on mixOmics (already installed in step 3)
# and a few additional packages listed below.

# pROC – needed for AUROC calculation inside MintTea
# (already installed in step 3 as pROC_1.18.5.tar.gz)
library(pROC)

# reshape2 – needed for MintTea output processing
# (already installed in step 3 as reshape2_1.4.4.tar.gz)
library(reshape2)

# mixOmics – core dependency of MintTea
library(mixOmics)


# ---- Install MintTea from GitHub ----
library(devtools)

cat("\nInstalling MintTea from GitHub: efratmuller/MintTea ...\n")
tryCatch({
  devtools::install_github(
    "efratmuller/MintTea",
    upgrade    = "never",
    build_vignettes = FALSE,
    quiet      = FALSE
  )
  library(MintTea)
  cat("\nMintTea installed and loaded successfully.\n")
  cat("MintTea version:", as.character(packageVersion("MintTea")), "\n")
}, error = function(e) {
  cat("\nERROR: Could not install MintTea from GitHub.\n")
  cat("Reason:", conditionMessage(e), "\n")
  cat("Please check your internet connection and try again.\n")
  cat("Manual installation command:\n")
  cat("  devtools::install_github('efratmuller/MintTea')\n")
})

print("MODULE MintTea INSTALLED")
