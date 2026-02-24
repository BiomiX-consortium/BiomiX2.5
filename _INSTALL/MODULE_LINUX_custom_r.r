# BiomiX – custom R packages not available on conda-forge/bioconda
# Covers: core packages, SNF-NEMO extension (NEMO), MintTea extension

script_dir <- dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE))
setwd(script_dir)

if (!file.exists("Package_linux")) {
  if (!file.exists("Package_linux.tar.gz")) {
    options(timeout = 6000)
    download.file(
      "https://github.com/IxI-97/BiomiX/releases/download/v2.4/Package_linux.tar.gz",
      destfile = "Package_linux.tar.gz", mode = "wb"
    )
  }
  untar("Package_linux.tar.gz")
}

setwd(file.path(script_dir, "Package_linux"))
options(Ncpus = 6)

# ── Core custom packages ───────────────────────────────────────────
remotes::install_local("masstools_1.0.13.tar.gz",  upgrade = "never")
remotes::install_local("massdataset_1.0.34.tar.gz", upgrade = "never")
remotes::install_local("metid_1.2.35.tar.gz",       upgrade = "never")
remotes::install_local("metpath_1.0.8.tar.gz",      upgrade = "never")
remotes::install_local("litsearchr_1.0.0.tar.gz",   upgrade = "never")
remotes::install_local("cmmr_1.0.5.tar.gz",         upgrade = "never")

# ── Methylation packages not available for R 4.4 on bioconda ─────
install.packages("ChAMPdata_2.36.0.tar.gz", repos = NULL, type = "source")
install.packages("ChAMP_2.34.0.tar.gz",     repos = NULL, type = "source")

# ── SNF/spatial package not on conda-forge ────────────────────────
install.packages("sabre_0.4.3.tar.gz", repos = NULL, type = "source")

# ── SNF-NEMO extension: NEMO ──────────────────────────────────────
install.packages("NEMO_0.1.0.tar.gz", repos = NULL, type = "source")

# ── MintTea extension (requires internet) ─────────────────────────
library(devtools)
tryCatch({
  devtools::install_github("efratmuller/MintTea", upgrade = "never",
                           build_vignettes = FALSE, quiet = FALSE)
  cat("MintTea installed successfully.\n")
}, error = function(e) {
  cat("WARNING: MintTea install failed (check internet connection).\n")
  cat("Reason:", conditionMessage(e), "\n")
})

library(masstools); library(massdataset); library(metid)
library(metpath);   library(litsearchr);  library(cmmr)
library(NEMO)

cat("Custom R packages installed successfully.\n")
