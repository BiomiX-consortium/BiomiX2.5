# BiomiX macOS – packages not available via conda on macOS
# Installs via BiocManager: minfi, MOFA2, and Illumina annotation packages
# that have unresolvable dependencies on macOS bioconda.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

options(Ncpus = 4, timeout = 600)

pkgs <- c(
  "minfi",
  "MOFA2",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kmanifest",
  "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  # Not on osx-arm64 bioconda
  "bsseq",
  "DSS",
  "DMRcate",
  "goseq"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    message("Already installed: ", pkg)
  }
}

message("MODULE_MAC_custom_r.r: done")
