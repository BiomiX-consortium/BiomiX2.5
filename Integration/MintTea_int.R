# MintTea multi-omics integration for BiomiX
# Reference: Muller et al., Nature Communications 2024
# https://doi.org/10.1038/s41467-024-46888-3
# Package: https://github.com/efratmuller/MintTea
#
# Identifies disease-associated multi-omic modules via repeated sparse gCCA
# (bootstrap-like consensus on top of DIABLO/mixOmics).
#
# # MANUAL INPUT (for testing outside BiomiX)
# args = list("CASE", "CTRL", "/path/to/BiomiX2.5")
# directory <- args[[3]]

# ---- Libraries ----

library(data.table)
library(vroom)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(rlist)
library(reshape2)
library(MintTea)
library(mixOmics)


# ---- Read configuration ----

MART <- vroom(paste(directory, "/Integration/x_BiomiX_DATABASE/mart_export_37.txt", sep = ""), delim = ",")
myList <- list()

COMMAND          <- vroom(paste(directory, "COMMANDS.tsv",          sep = "/"), delim = "\t")
COMMAND_MOFA     <- vroom(paste(directory, "COMMANDS_MOFA.tsv",     sep = "/"), delim = "\t")
COMMAND_ADVANCED <- vroom(paste(directory, "COMMANDS_ADVANCED.tsv", sep = "/"), delim = "\t")

Max_features <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[3])

# ---- MintTea parameters ----

n_repeats      <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MINTTEA_OPTIONS[1])
n_folds        <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MINTTEA_OPTIONS[2])
edge_threshold <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MINTTEA_OPTIONS[3])
keepX_per_comp <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MINTTEA_FEATURES[1])
design_val     <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MINTTEA_FEATURES[2])

DIR_METADATA <- readLines(paste(directory, "directory.txt", sep = "/"))

if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
  METADATA <- read_excel(DIR_METADATA)
  print("Metadata Excel File read successfully!")
} else {
  METADATA <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)
}


# ---- Data processing functions (shared with DIABLO/SNF pipelines) ----

Undefined_processing <- function(matrix, mer) {
  matrix <- matrix %>% filter(CONDITION == args[2] | CONDITION == args[1])
  lst    <- list()
  sample <- matrix$ID
  matrix <- matrix %>% dplyr::select(!c(CONDITION, ID))
  lst[[1]] <- as.data.frame(matrix)
  rownames(lst[[1]]) <- sample
  lst[[2]] <- colnames(matrix)
  return(lst)
}

Metabolomics_processing <- function(annotation, matrix, mer) {
  annotation <- annotation %>% distinct(NAME.x, .keep_all = TRUE)
  annotation$NAMES <- paste(annotation$NAME.x, annotation$Name, sep = "/")
  x      <- as.data.frame(colnames(matrix))
  x$id   <- 1:nrow(x)
  colnames(x)[1] <- "Peaks"
  xx <- merge(x, annotation, by.x = "Peaks", by.y = "NAME.x", all.x = TRUE)
  xx <- xx[order(xx$id), ]
  for (i in 1:nrow(xx)) {
    if (!is.na(xx$NAMES[i])) xx$Peaks[i] <- xx$NAMES[i]
  }
  serum_metabolomics_EASY <- xx$Peaks
  lst    <- list(annotation)
  matrix <- matrix %>% filter(CONDITION == args[2] | CONDITION == args[1])
  sample <- matrix$ID
  matrix <- matrix %>% dplyr::select(!c(CONDITION, ID))
  lst[[1]] <- as.data.frame(matrix)
  rownames(lst[[1]]) <- sample
  lst[[2]] <- serum_metabolomics_EASY
  return(lst)
}

Transcriptomics_processing <- function(annotation, matrix, mer) {
  lst    <- list()
  matrix <- merge(matrix, MART, by.x = "ID", by.y = "Gene.stable.ID")
  matrix <- matrix %>% arrange(desc(variance))
  matrix <- matrix[1:Max_features, ]
  Bcell_RNAseq_VIEW <- as.character(matrix$`Gene name`)
  matrix <- t(matrix)
  colnames(matrix) <- matrix[1, ]
  matrix <- as.data.frame(matrix)
  matrix <- matrix[c(-1, -(nrow(matrix) - 1), -nrow(matrix)), ]
  Bcell_RNAseq_VIEW <- make.unique(Bcell_RNAseq_VIEW, sep = "_")
  colnames(matrix) <- Bcell_RNAseq_VIEW
  lst[[1]] <- matrix
  lst[[2]] <- Bcell_RNAseq_VIEW
  return(lst)
}

Methylomics_processing <- function(annotation, matrix, metadata, mer) {
  annotation <- annotation[, c("gene", "CpG_island")]
  matrix     <- merge(matrix, annotation, by.x = "ID", by.y = "CpG_island")
  matrix     <- matrix[1:Max_features, ]
  Methylome_WB_VIEW <- as.character(matrix$gene)
  Methylome_WB_EASY <- paste(matrix$ID, matrix$gene, sep = "/")
  lst    <- list()
  matrix <- t(matrix)
  colnames(matrix) <- matrix["ID", ]
  matrix <- as.data.frame(matrix)
  y      <- rownames(matrix) %in% c("ID", "variance")
  matrix <- matrix[!y, ]
  matrix <- matrix[c(-(nrow(matrix) - 1), -nrow(matrix)), ]
  Methylome_WB_EASY <- make.unique(Methylome_WB_EASY, sep = "_")
  colnames(matrix)  <- Methylome_WB_EASY
  lst[[1]] <- matrix
  lst[[2]] <- Methylome_WB_EASY
  return(lst)
}


# ---- Load input data ----

myList   <- list()
names_X  <- c()

for (i in 1:length(COMMAND$INTEGRATION)) {

  if (COMMAND$INTEGRATION[i] == "YES") {

    if (COMMAND$DATA_TYPE[i] == "Metabolomics") {
      directory2      <- paste(directory, "/Integration/INPUT/Metabolomics_", COMMAND$LABEL[i], "_", args[1], "_vs_", args[2], sep = "")
      serum_metab     <- vroom(paste(directory2, "/Metabolomics_", COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
      directory2      <- paste(directory, "/Metabolomics/OUTPUT/", COMMAND$LABEL[i], "_", args[1], "_vs_", args[2], sep = "")
      serum_annot     <- vroom(paste(directory2, "/", COMMAND$LABEL[i], "_", args[1], "_vs_", args[2], "_results.tsv", sep = ""), delim = "\t")
      INPUTX          <- Metabolomics_processing(serum_annot, serum_metab, COMMAND$LABEL[i])
      myList          <- list.append(myList, INPUTX[[1]])
      names_X         <- append(names_X, COMMAND$LABEL[i])
    }

    if (COMMAND$DATA_TYPE[i] == "Transcriptomics") {
      directory2      <- paste(directory, "/Integration/INPUT/", COMMAND$LABEL[i], "_", args[1], "_vs_", args[2], sep = "")
      rna_matrix      <- vroom(paste(directory2, "/", COMMAND$LABEL[i], "_", args[1], "_vs_", args[2], "_normalized_vst_variance.tsv", sep = ""), delim = "\t")
      rna_meta        <- vroom(paste(directory2, "//Metadata_", COMMAND$LABEL[i], "_", args[1], ".tsv", sep = ""), delim = "\t")
      INPUTX          <- Transcriptomics_processing(rna_meta, rna_matrix, COMMAND$LABEL[i])
      myList          <- list.append(myList, INPUTX[[1]])
      names_X         <- append(names_X, COMMAND$LABEL[i])
    }

    if (COMMAND$DATA_TYPE[i] == "Methylomics") {
      directory2      <- paste(directory, "/Integration/INPUT/Methylome_", COMMAND$LABEL[i], "_", args[1], "_vs_", args[2], sep = "")
      methy_matrix    <- vroom(paste(directory2, "/", COMMAND$LABEL[i], "_matrix_MOFA.tsv", sep = ""), delim = "\t")
      methy_meta      <- vroom(paste(directory2, "/", COMMAND$LABEL[i], "_metadata_MOFA.tsv", sep = ""), delim = "\t")
      directory2      <- paste(directory, "/Methylomics/OUTPUT/", COMMAND$LABEL[i], "_", args[1], "_vs_", args[2], sep = "")
      methy_annot     <- vroom(paste(directory2, "/DMP_", COMMAND$LABEL[i], "_Methylome_", args[1], "_vs_", args[2], ".tsv", sep = ""), delim = "\t", col_names = TRUE)
      INPUTX          <- Methylomics_processing(methy_annot, methy_matrix, methy_meta, COMMAND$LABEL[i])
      myList          <- list.append(myList, INPUTX[[1]])
      names_X         <- append(names_X, COMMAND$LABEL[i])
    }

    if (COMMAND$DATA_TYPE[i] == "Undefined") {
      directory2      <- paste(directory, "/Integration/INPUT/Undefined_", COMMAND$LABEL[i], "_", args[1], "_vs_", args[2], sep = "")
      undef_matrix    <- vroom(paste(directory2, "/Undefined_", COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
      INPUTX          <- Undefined_processing(undef_matrix, COMMAND$LABEL[i])
      myList          <- list.append(myList, INPUTX[[1]])
      names_X         <- append(names_X, COMMAND$LABEL[i])
    }
  }
}

names(myList) <- names_X


# ---- Find common samples across all omics ----

common_rows <- Reduce(function(x, y) intersect(rownames(x), rownames(y)), myList)
cat("Common samples across all omics:", length(common_rows), "\n")


# ---- Build MintTea input: single data.frame with omic-prefixed columns ----
# MintTea requires columns named as "<prefix>__<feature>", e.g. "RNA__TP53"
# and a "label" column with "case" / "control" values.

prefixed_list <- lapply(seq_along(myList), function(i) {
  df  <- as.data.frame(myList[[i]])
  df  <- df[common_rows, , drop = FALSE]
  mat <- apply(df, 2, as.numeric)
  df  <- as.data.frame(mat)
  rownames(df) <- common_rows
  # Remove CV-fold-unsafe features: a feature is risky if the number of
  # samples that differ from the modal value is ≤ ceil(N/n_folds).
  # Such features can become zero-variance in a CV training fold and
  # cause block.splsda to error.  ceil(N/n_folds) is the largest held-out
  # fold size; requiring minority_count > that guarantees safety.
  n_max_holdout  <- ceiling(nrow(df) / n_folds)
  feat_modal     <- apply(df, 2, function(x) max(table(x)))
  feat_minority  <- nrow(df) - feat_modal
  df <- df[, feat_minority > n_max_holdout, drop = FALSE]
  # Prefix column names with omic label, sanitise special characters
  colnames(df) <- paste0(names_X[i], "__", make.names(colnames(df)))
  df
})

combined_data            <- as.data.frame(do.call(cbind, prefixed_list))
rownames(combined_data)  <- common_rows

# Align metadata and add required MintTea columns:
#   sample_id     – sample identifier
#   disease_state – "disease" (case) or "healthy" (control)
METADATA_f <- METADATA[METADATA$ID %in% common_rows, ]
METADATA_f <- METADATA_f[order(factor(METADATA_f$ID, levels = common_rows)), ]
combined_data$sample_id     <- common_rows
combined_data$disease_state <- ifelse(METADATA_f$CONDITION == args[1], "disease", "healthy")

cat("MintTea input dimensions:", nrow(combined_data), "samples x",
    ncol(combined_data) - 2, "features\n")
cat("Label distribution:\n")
print(table(combined_data$disease_state))


# ---- Create output directory ----

out_dir <- paste(directory, "/Integration/OUTPUT/MintTea_", args[1], "_vs_", args[2], sep = "")
dir.create(path = out_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
setwd(out_dir)


# ---- Run MintTea ----

cat("\n\n--- Running MintTea ---\n")
cat("View prefixes     :", paste(names_X, collapse = ", "), "\n")
cat("n_repeats         :", n_repeats, "\n")
cat("n_folds           :", n_folds, "\n")
cat("edge_threshold    :", edge_threshold, "\n")
cat("keepX per comp    :", keepX_per_comp, "\n")
cat("design value      :", design_val, "\n\n")

minttea_results <- MintTea(
  proc_data             = combined_data,
  view_prefixes         = names_X,
  sample_id_column      = "sample_id",
  study_group_column    = "disease_state",
  case_group_name       = "disease",
  control_group_name    = "healthy",
  param_n_repeats       = c(n_repeats),
  param_n_folds         = c(n_folds),
  param_edge_thresholds = c(edge_threshold),
  param_diablo_keepX    = c(keepX_per_comp),
  param_sgcca_design    = c(design_val)
)

# Save raw results object
saveRDS(minttea_results,
        file = paste0("MintTea_", args[1], "_vs_", args[2], "_results.rds"))

n_modules <- length(minttea_results)
cat("\n--- MintTea finished:", n_modules, "module(s) identified ---\n")


# ---- Downstream analysis and output generation ----

if (n_modules == 0) {

  cat("No modules were identified.\n")
  cat("Try lowering edge_threshold, increasing n_repeats, or lowering keepX.\n")

} else {

  # 1. Module summary table ----
  summary_df <- do.call(rbind, lapply(seq_len(n_modules), function(m) {
    mod <- minttea_results[[m]]
    data.frame(
      module                   = m,
      module_size              = mod$module_size,
      auroc                    = round(mod$auroc, 4),
      shuffled_auroc_mean      = round(mean(mod$shuffled_auroc), 4),
      inter_view_corr          = round(mod$inter_view_corr, 4),
      shuffled_inter_view_corr = round(mean(mod$shuffled_inter_view_corr), 4)
    )
  }))

  write.table(summary_df,
    file      = paste0("MintTea_", args[1], "_vs_", args[2], "_module_summary.tsv"),
    quote     = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  cat("Module summary saved.\n")


  # 2. Per-module feature tables ----
  for (m in seq_len(n_modules)) {
    mod <- minttea_results[[m]]
    for (view in names(mod$features)) {
      feats <- mod$features[[view]]
      if (length(feats) > 0) {
        feat_df <- data.frame(feature = feats, view = view, module = m)
        write.table(feat_df,
          file  = paste0("MintTea_module", m, "_", view, "_features.tsv"),
          quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
      }
    }
  }
  cat("Feature tables saved.\n")


  # 3. Edge (co-occurrence) tables ----
  for (m in seq_len(n_modules)) {
    mod <- minttea_results[[m]]
    if (!is.null(mod$module_edges) && nrow(mod$module_edges) > 0) {
      write.table(mod$module_edges,
        file  = paste0("MintTea_module", m, "_edges.tsv"),
        quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }
  }
  cat("Edge tables saved.\n")


  # 4. PDF plots ----
  pdf(file   = paste0("MintTea_", args[1], "_vs_", args[2], "_plots.pdf"),
      width  = 12, height = 8)


  # 4a. AUROC bar chart: observed vs shuffled baseline ----
  auroc_data <- do.call(rbind, lapply(seq_len(n_modules), function(m) {
    mod <- minttea_results[[m]]
    rbind(
      data.frame(module = paste0("Module ", m), value = mod$auroc,
                 type   = "Observed"),
      data.frame(module = paste0("Module ", m),
                 value  = mean(mod$shuffled_auroc),
                 type   = "Shuffled (mean)")
    )
  }))

  p_auroc <- ggplot(auroc_data, aes(x = module, y = value, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Observed" = "steelblue",
                                 "Shuffled (mean)" = "grey70")) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
    labs(title  = paste("MintTea Module AUROC –", args[1], "vs", args[2]),
         x      = "Module", y = "AUROC", fill = "") +
    theme_classic() +
    theme(axis.text       = element_text(size = 11),
          plot.title      = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom")
  print(p_auroc)


  # 4b. Inter-view correlation bar chart ----
  corr_data <- do.call(rbind, lapply(seq_len(n_modules), function(m) {
    mod <- minttea_results[[m]]
    rbind(
      data.frame(module = paste0("Module ", m), value = mod$inter_view_corr,
                 type   = "Observed"),
      data.frame(module = paste0("Module ", m),
                 value  = mean(mod$shuffled_inter_view_corr),
                 type   = "Shuffled (mean)")
    )
  }))

  p_corr <- ggplot(corr_data, aes(x = module, y = value, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Observed" = "darkorange",
                                 "Shuffled (mean)" = "grey70")) +
    labs(title  = paste("MintTea Inter-view Correlation –", args[1], "vs", args[2]),
         x      = "Module", y = "Inter-view Correlation", fill = "") +
    theme_classic() +
    theme(axis.text       = element_text(size = 11),
          plot.title      = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom")
  print(p_corr)


  # 4c. Module overview scatter plot (AUROC vs inter-view corr, size = n features) ----
  p_overview <- ggplot(summary_df,
                       aes(x = inter_view_corr, y = auroc, size = module_size,
                           label = paste0("M", module))) +
    geom_point(colour = "steelblue", alpha = 0.7) +
    geom_text(nudge_y = 0.015, size = 3.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red", alpha = 0.5) +
    labs(title  = paste("MintTea Module Overview –", args[1], "vs", args[2]),
         x      = "Inter-view Correlation",
         y      = "AUROC",
         size   = "Module Size") +
    theme_classic() +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  print(p_overview)


  # 4d. Per-module feature bar charts by omic ----
  for (m in seq_len(n_modules)) {
    mod  <- minttea_results[[m]]
    all_feats <- do.call(rbind, lapply(names(mod$features), function(view) {
      feats <- mod$features[[view]]
      if (length(feats) > 0) data.frame(feature = feats, view = view)
      else NULL
    }))
    if (!is.null(all_feats) && nrow(all_feats) > 0) {
      p_feat <- ggplot(all_feats,
                       aes(x = reorder(feature, as.numeric(as.factor(view))),
                           fill = view)) +
        geom_bar(stat = "count") +
        coord_flip() +
        labs(title  = paste("Module", m, "features  |  AUROC:", round(mod$auroc, 3)),
             x      = "Feature", y = "Count", fill = "Omic") +
        theme_classic() +
        theme(axis.text.y = element_text(size = 8))
      print(p_feat)
    }
  }


  # 4e. Feature co-occurrence heatmap per module ----
  for (m in seq_len(n_modules)) {
    mod   <- minttea_results[[m]]
    edges <- mod$module_edges
    if (!is.null(edges) && nrow(edges) > 0 &&
        all(c("feature1", "feature2", "weight") %in% colnames(edges))) {

      # Limit to top 30 feature pairs for readability
      edges <- edges[order(-edges$weight), ]
      top_feats <- unique(c(as.character(edges$feature1[1:min(30, nrow(edges))]),
                            as.character(edges$feature2[1:min(30, nrow(edges))])))
      edges_sub <- edges[edges$feature1 %in% top_feats &
                         edges$feature2 %in% top_feats, ]

      all_feat <- unique(c(as.character(edges_sub$feature1),
                           as.character(edges_sub$feature2)))
      mat <- matrix(0, nrow = length(all_feat), ncol = length(all_feat),
                    dimnames = list(all_feat, all_feat))
      for (r in seq_len(nrow(edges_sub))) {
        f1 <- as.character(edges_sub$feature1[r])
        f2 <- as.character(edges_sub$feature2[r])
        w  <- edges_sub$weight[r]
        mat[f1, f2] <- w
        mat[f2, f1] <- w
      }

      mat_long <- melt(mat)
      p_heat <- ggplot(mat_long, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(title = paste("Module", m, "– feature co-occurrence"),
             fill  = "Weight") +
        theme_classic() +
        theme(axis.text.x  = element_text(angle = 90, hjust = 1, size = 6),
              axis.text.y  = element_text(size = 6),
              axis.title   = element_blank())
      print(p_heat)
    }
  }

  dev.off()
  cat("PDF plots saved.\n")

  cat("\n--- MintTea output written to:", out_dir, "---\n")
}
