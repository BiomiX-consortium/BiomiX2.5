# ============================================================
# MintTea real-data test using the MTBLS7623 PTB dataset
# (Transcriptomics + Metabolomics, PTB vs Healthy Controls)
#
# Run this script from R or RStudio on any machine.
# No BiomiX installation required.
#
# Reference: Muller et al., Nat Commun 2024
# https://github.com/efratmuller/MintTea
# ============================================================

library(MintTea)
library(vroom)
library(dplyr)
library(ggplot2)

# ---- Paths (edit if your BiomiX2.5 folder is in a different location) ----

# Auto-detect script location (works in radian, base R console, and Rscript)
base_dir <- tryCatch({
  # When sourced in radian / base R console
  dirname(normalizePath(sys.frame(1)$ofile))
}, error = function(e) tryCatch({
  # When run via: Rscript test_MintTea_real_data.R
  args <- commandArgs(trailingOnly = FALSE)
  script_flag <- args[grep("--file=", args)]
  if (length(script_flag) > 0)
    dirname(normalizePath(sub("--file=", "", script_flag)))
  else
    getwd()
}, error = function(e) getwd()))

# If auto-detection fails, uncomment and set manually:
# base_dir <- "/Users/lourdes/Dropbox/Curro/BiomiX2.5"

cat("Base directory:", base_dir, "\n")

rna_file   <- file.path(base_dir, "Transcriptomics/INPUT/MTBLS7623_RNA_seq.tsv")
metab_file <- file.path(base_dir, "Metabolomics/INPUT/MTBLS7623_positive_mode_unannotated.tsv")
meta_file  <- file.path(base_dir, "Example_dataset/MTBLS7623-PRJNA971365/Metadata/MTBLS7623_Metadata.tsv")


# ---- Parameters ----

N_GENES    <- 100    # top variable genes to keep (reduce to speed up)
CASE       <- "PTB"  # case group (also try "PTB_DM")
CONTROL    <- "HC"   # control group


# ---- 1. Load metadata and filter to binary comparison ----

cat("Loading metadata...\n")
metadata <- vroom(meta_file, delim = "\t", col_select = c("ID", "CONDITION"))
metadata <- metadata %>% filter(CONDITION %in% c(CASE, CONTROL))
samples  <- metadata$ID
cat("Samples:", length(samples), "(", CASE, ":", sum(metadata$CONDITION == CASE),
    "| ", CONTROL, ":", sum(metadata$CONDITION == CONTROL), ")\n")


# ---- 2. Load and prepare RNA-seq ----
# File format: genes x samples (needs transposing)
# Values are raw counts → log1p-normalise, then keep top variable genes

cat("Loading RNA-seq...\n")
rna_raw <- vroom(rna_file, delim = "\t")
gene_ids <- rna_raw$ID
rna_raw  <- rna_raw %>% dplyr::select(-ID)

# Keep only the samples in our comparison
rna_raw  <- rna_raw %>% dplyr::select(any_of(samples))

# Transpose: rows = samples, cols = genes
rna_mat  <- as.data.frame(t(rna_raw))
colnames(rna_mat) <- gene_ids
rna_mat  <- rna_mat[samples[samples %in% rownames(rna_mat)], ]

# log1p normalise
rna_mat  <- log1p(apply(rna_mat, 2, as.numeric))
rna_mat  <- as.data.frame(rna_mat)
rownames(rna_mat) <- samples[samples %in% rownames(as.data.frame(t(rna_raw)))]

# Select top N_GENES by variance
gene_var <- apply(rna_mat, 2, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:N_GENES]
rna_mat  <- rna_mat[, top_genes, drop = FALSE]
cat("RNA matrix:", nrow(rna_mat), "samples x", ncol(rna_mat), "genes\n")


# ---- 3. Load and prepare metabolomics ----
# File format: peaks x samples (needs transposing)
# Values are intensities → log1p-normalise

cat("Loading metabolomics...\n")
metab_raw  <- vroom(metab_file, delim = "\t")
peak_ids   <- metab_raw$ID
metab_raw  <- metab_raw %>% dplyr::select(-ID)

# Keep only the samples in our comparison
metab_raw  <- metab_raw %>% dplyr::select(any_of(samples))

# Transpose: rows = samples, cols = peaks
metab_mat  <- as.data.frame(t(metab_raw))
colnames(metab_mat) <- peak_ids

# log1p normalise
metab_mat  <- log1p(apply(metab_mat, 2, as.numeric))
metab_mat  <- as.data.frame(metab_mat)
metab_samples <- samples[samples %in% rownames(as.data.frame(t(metab_raw)))]
rownames(metab_mat) <- metab_samples
cat("Metabolomics matrix:", nrow(metab_mat), "samples x", ncol(metab_mat), "peaks\n")


# ---- 4. Align samples across both omics ----

common_samples <- intersect(rownames(rna_mat), rownames(metab_mat))
rna_mat   <- rna_mat[common_samples, , drop = FALSE]
metab_mat <- metab_mat[common_samples, , drop = FALSE]
metadata  <- metadata[metadata$ID %in% common_samples, ]
metadata  <- metadata[order(match(metadata$ID, common_samples)), ]
cat("Common samples:", length(common_samples), "\n")


# ---- 4b. Remove CV-fold-unsafe features ----
# A feature is "safe" only if enough samples deviate from the modal value to
# guarantee non-zero variance in every CV training fold. With vfold_cv(v=5)
# over N samples, the largest test fold holds ceil(N/5) samples.  If ALL
# minority-value samples fit in that one held-out fold, the training set is
# left with a constant feature (zero variance) → block.splsda error.
# Requiring minority_count > ceil(N/5) guarantees safety.

n_max_holdout <- ceiling(length(common_samples) / 5)  # e.g. ceil(26/5) = 6

cv_safe <- function(mat) {
  apply(mat, 2, function(x) {
    modal_count    <- max(table(x))
    minority_count <- length(x) - modal_count
    minority_count > n_max_holdout
  })
}

rna_mat  <- rna_mat[,   cv_safe(rna_mat),   drop = FALSE]
cat("RNA features after CV-safe filter:", ncol(rna_mat), "\n")

metab_mat <- metab_mat[, cv_safe(metab_mat), drop = FALSE]
cat("Metabolomics features after CV-safe filter:", ncol(metab_mat), "\n")


# ---- 5. Build MintTea input ----
# Requires a single data.frame with:
#   - omic-prefixed columns: "RNA__gene", "METAB__peak"
#   - a sample_id column
#   - a disease_state column ("disease" / "healthy")

colnames(rna_mat)   <- paste0("RNA__",   make.names(colnames(rna_mat)))
colnames(metab_mat) <- paste0("METAB__", make.names(colnames(metab_mat)))

combined <- cbind(rna_mat, metab_mat)
combined$sample_id     <- common_samples
combined$disease_state <- ifelse(metadata$CONDITION == CASE, "disease", "healthy")

cat("\nLabel distribution:\n")
print(table(combined$disease_state))
cat("Total features:", ncol(combined) - 2, "\n")


# ---- 6. Run MintTea ----

cat("\n--- Running MintTea (", CASE, "vs", CONTROL, ") ---\n")

set.seed(42)
results <- MintTea(
  proc_data             = combined,
  view_prefixes         = c("RNA", "METAB"),
  sample_id_column      = "sample_id",
  study_group_column    = "disease_state",
  case_group_name       = "disease",
  control_group_name    = "healthy",
  param_n_repeats       = c(20),
  param_n_folds         = c(5),
  param_edge_thresholds = c(0.5),
  param_diablo_keepX    = c(10),
  param_sgcca_design    = c(0.5)
)

# results structure: results[[run_id]][[module_id]] = module_data
# (MintTea returns a 2-level nested list; the top level is the parameter
#  setting string, the second level is "module1", "module2", etc.)

if (length(results) == 0) {
  cat("\n--- MintTea finished: 0 module(s) found ---\n")
  cat("No modules identified. Try lowering param_edge_thresholds (e.g. 0.7)\n")
  cat("or increasing param_n_repeats.\n")
} else {
  run_modules <- results[[1]]          # first (and only) parameter setting
  n_modules   <- length(run_modules)
  cat("\n--- MintTea finished:", n_modules, "module(s) found ---\n")


  # ---- 7. Print results summary ----

  for (m in seq_len(n_modules)) {
    mod <- run_modules[[m]]
    cat("\n=== Module", m, "===\n")
    cat("  Size (features)  :", mod$module_size, "\n")
    cat("  AUROC            :", round(mod$auroc, 3),
        "(shuffled mean:", round(mean(mod$shuffled_auroc), 3), ")\n")
    cat("  Inter-view corr  :", round(mod$inter_view_corr, 3),
        "(shuffled mean:", round(mean(mod$shuffled_inter_view_corr), 3), ")\n")
    # mod$features is a flat character vector like c("RNA__GENE1", "METAB__peak1")
    cat("  Features by omic :\n")
    feats <- mod$features
    for (pfx in c("RNA", "METAB")) {
      view_feats <- grep(paste0("^", pfx, "__"), feats, value = TRUE)
      if (length(view_feats) > 0) {
        clean <- sub(paste0("^", pfx, "__"), "", view_feats)
        cat("   ", pfx, ":", paste(head(clean, 5), collapse = ", "),
            if (length(clean) > 5) paste0("... (+", length(clean) - 5, " more)") else "", "\n")
      }
    }
  }

  # ---- 8. Quick AUROC plot ----
  auroc_df <- do.call(rbind, lapply(seq_len(n_modules), function(m) {
    mod <- run_modules[[m]]
    rbind(
      data.frame(module = paste0("M", m), value = mod$auroc,               type = "Observed"),
      data.frame(module = paste0("M", m), value = mean(mod$shuffled_auroc), type = "Shuffled (mean)")
    )
  }))

  p <- ggplot(auroc_df, aes(x = module, y = value, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Observed" = "steelblue", "Shuffled (mean)" = "grey70")) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
    labs(title  = paste("MintTea –", CASE, "vs", CONTROL),
         x = "Module", y = "AUROC", fill = "") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}
