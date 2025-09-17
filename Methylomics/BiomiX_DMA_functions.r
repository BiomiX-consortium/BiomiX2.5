# Functions for the methylomics pipeline BiomiX_DMA.r script



filter_metadata_matrix <- function(COMMAND_ADVANCED, Metadata, Matrix) {
  # --- Code from your original script (unchanged) ---
  METADATA_FILT <- !is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])
  METADATA_FILT_INDEX <-grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))
  
  repetition = 0
  for (meta_filter in METADATA_FILT_INDEX){
    repetition <- repetition + 1 
    if (!is.na(COMMAND_ADVANCED[3,grep( "*.FILTERING.*", colnames(COMMAND_ADVANCED))])[repetition]){
      COLNAME<-as.character(COMMAND_ADVANCED[1,meta_filter])
      if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Numerical"){
        To_filter<-as.numeric(unlist(Metadata[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,1)
        characters_to_remove <- c(">", "<", "=", " ")
        value_threshold <- as.numeric(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "<" = function(a, b) a < b,
                                      ">" = function(a, b) a > b,
                                      "=" = function(a, b) a == b,
                                      ">=" = function(a, b) a >= b,
                                      "<=" = function(a, b) a <= b,
                                      NA)
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold),]
        Matrix <- Matrix[,comparison_operator(To_filter, value_threshold)]
        
      }else if (as.character(COMMAND_ADVANCED[2,meta_filter]) =="Factors"){
        To_filter<- as.character(unlist(Metadata[,COLNAME]))
        simbol<-substr(as.character(COMMAND_ADVANCED[3,meta_filter]),1,2)
        characters_to_remove <- c("!=", "==", " ")
        value_threshold <- as.character(gsub(paste(characters_to_remove, collapse = "|"), "", as.character(COMMAND_ADVANCED[3,meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "==" = function(a, b) a == b,
                                      "!=" = function(a, b) a != b,
                                      NA)
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold),]
        Matrix <- Matrix[,comparison_operator(To_filter, value_threshold)]
      }
    }
  }
  
  # --- Make sure the variables exist outside the function (assign to global env) ---
  assign("METADATA_FILT", METADATA_FILT, envir = .GlobalEnv)
  assign("METADATA_FILT_INDEX", METADATA_FILT_INDEX, envir = .GlobalEnv)
  assign("Metadata", Metadata, envir = .GlobalEnv)
  assign("Matrix", Matrix, envir = .GlobalEnv)

}


safe_transpose_df <- function(df) {
  # Step 1: Ensure all columns are numeric before transpose
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(as.character(x))))
  
  # Step 2: Transpose and convert back to data.frame
  df_t <- as.data.frame(t(as.matrix(df)), stringsAsFactors = FALSE)
  
  # Step 3: Convert back to numeric after transpose
  df_t[] <- lapply(df_t, function(x) suppressWarnings(as.numeric(as.character(x))))
  
  return(df_t)
}



qc_preview_visualization <- function(Matrix, Metadata_individual, directory) {
  # --- Original code preserved ---
  Samples_preview <- colnames(Matrix[,-1])
  cpg <- Matrix$ID
  numeric_data <- Matrix[,-1]
  numeric_data <- transpose(numeric_data)
  numeric_data <- apply(numeric_data, 2, as.numeric)
  rownames(numeric_data) <- Samples_preview
  colnames(numeric_data) <- cpg
  
  metadata <- Metadata_individual
  numeric_data <- as.data.frame(numeric_data)
  colnames(metadata)[colnames(metadata) == "CONDITION"] <- "condition"
  
  # Load the BiomiX preview Shiny app
  source(paste(directory,'/BiomiX_preview.r', sep=""))
  
  browser_analysis <- readLines(paste(directory,'/_INSTALL/CHOISE_BROWSER_pre-view', sep=""), n = 1)
  
  options(browser = get_default_browser())
  print("Pre QC data visualization. If the browser does not open automatically copy and paste the link on a browser")
  print("Once completed, click on close app to continue the analysis")
  print("ATTENTION!: If the browser does not open, set the path to your browser in the CHOISE_BROWSER_pre-view file in the BiomiX _INSTALL folder")
  
  Preview <- shiny::runApp(runShinyApp(numeric_data, metadata), launch.browser = TRUE)
  
  Matrix <- as.data.frame(t(Preview$matrix))
  Metadata_individual <- Preview$metadata
  colnames(Metadata_individual)[colnames(Metadata_individual) == "condition"] <- "CONDITION"
  
  Matrix <- cbind(cpg, Matrix)
  colnames(Matrix)[1] <- "ID"
  
  # --- Return updated objects ---
  assign("Matrix", Matrix, envir = .GlobalEnv)
  assign("Metadata_individual", Metadata_individual, envir = .GlobalEnv)
  assign("Samples_preview", Samples_preview, envir = .GlobalEnv)
  assign("cpg", cpg, envir = .GlobalEnv)
  assign("numeric_data", numeric_data, envir = .GlobalEnv)
  assign("metadata", metadata, envir = .GlobalEnv)
  assign("browser_analysis", browser_analysis, envir = .GlobalEnv)
  assign("Preview", Preview, envir = .GlobalEnv)
}




perform_statistical_dmp_analysis <- function(Matrix,
                                             Metadata_individual,
                                             directory,
                                             Cell_type,
                                             args,
                                             array,
                                             padju,
                                             LogFC) {
  # Create output directories
  dir.create(path = file.path(directory, "Methylomics/OUTPUT/"), showWarnings = FALSE)
  dir.create(path = file.path(directory, "Methylomics/OUTPUT/", paste0(Cell_type, "_", args[1], "_vs_", args[2])), 
             showWarnings = FALSE)
  
  directory2 <- file.path(directory, "Methylomics/OUTPUT/", paste0(Cell_type, "_", args[1], "_vs_", args[2]))
  setwd(directory2)
  
  library(vroom)
  library(dplyr)
  library(ChAMP)
  
  spot <- Matrix$ID
  Matrix <- as.data.frame(Matrix[, -1])
  rownames(Matrix) <- spot
  
  Metadata_individual$CONDITION <- as.factor(Metadata_individual$CONDITION)
  
  Matrix <- as.matrix(Matrix)
  Matrix <- apply(Matrix, 2, as.numeric)
  Matrix <- as.data.frame(Matrix)
  rownames(Matrix) <- spot
  
  myDMP <- champ.DMP(beta = Matrix, pheno = Metadata_individual$CONDITION, adjPVal = 1, arraytype = array)
  
  results <- myDMP[[1]]
  results$CpG_island <- rownames(results)
  
  setwd(directory2)
  
  write.table(results, 
              file = paste0("DMP_", Cell_type, "_Methylome_", args[1], "_vs_", args[2], ".tsv"),  
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  ALTI <- subset(results, adj.P.Val < padju)
  ALTI <- subset(ALTI, deltaBeta > LogFC)
  ALTI <- ALTI[order(ALTI$adj.P.Val), ]
  
  BASSI <- subset(results, adj.P.Val < padju)
  BASSI <- subset(BASSI, deltaBeta < -LogFC)
  BASSI <- BASSI[order(BASSI$adj.P.Val), ]
  
  NO <- subset(results, deltaBeta < LogFC & deltaBeta > -LogFC)
  
  pdf(file = paste0("plot_DMA_", args[1], "_", args[2], "_", Cell_type, ".pdf"))
  
  library(ggplot2)
  library(ggrepel)
  
  p <- ggplot() +
    ggtitle(paste0("DMA_", Cell_type, "_", args[1], "_vs_", args[2])) + 
    theme(
      plot.title = element_text(color = "black", size = 14, face = "bold.italic", hjust = 0.5),
      axis.title.x = element_text(color = "black", size = 14, face = "bold"),
      axis.title.y = element_text(color = "black", size = 14, face = "bold")
    ) +
    geom_point(data = ALTI, aes(x = deltaBeta, y = -log10(adj.P.Val)), color = "red") +
    geom_text_repel(data = ALTI[c(1:25), ], aes(x = deltaBeta , y = -log10(adj.P.Val), label=gene),
                    hjust = 0, vjust = 0, size = 3, arrow = NULL, force = 10, force_pull = 1, max.overlaps = 25) +
    geom_point(data = BASSI, aes(x = deltaBeta, y = -log10(adj.P.Val)), color = "blue") +
    geom_text_repel(data = BASSI[1:25, ], aes(x = deltaBeta , y = -log10(adj.P.Val), label=gene),
                    hjust = 0, vjust = 0, size = 3, force = 20, force_pull = 1, max.overlaps = 25) +
    geom_point(data = NO, aes(x = deltaBeta, y = -log10(adj.P.Val)), color = "black")
  
  p2 <- p + 
    geom_vline(xintercept = c(-LogFC, LogFC), col = "red") +
    geom_hline(yintercept = -log10(padju), col = "red")
  
  print(p2)
  dev.off()
  
  write.table(ALTI$gene, file = paste0(Cell_type,"_",args[1],"_",args[2],"_UP.tsv"), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(BASSI$gene, file = paste0(Cell_type,"_",args[1],"_",args[2],"_DOWN.tsv"), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(ALTI, file = paste0(Cell_type,"_",args[1],"_",args[2],"_GENESUP.tsv"), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(BASSI, file = paste0(Cell_type,"_",args[1],"_",args[2],"_GENESDOWN.tsv"), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(Metadata_individual[order(Metadata_individual$CONDITION, decreasing = TRUE),],
              file = paste0(args[1], "_", args[2], "_", Cell_type, "_metadata.tsv"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  dir.create(path = paste0(directory2, "/Pathway_analysis"), 
             showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  directory_path <- paste0(directory2, "/Pathway_analysis")
  
  # Export to global environment to keep compatibility with downstream code
  assign("directory2", directory2, envir = .GlobalEnv)
  assign("Matrix", Matrix, envir = .GlobalEnv)
  assign("results", results, envir = .GlobalEnv)
  assign("ALTI", ALTI, envir = .GlobalEnv)
  assign("BASSI", BASSI, envir = .GlobalEnv)
  assign("NO", NO, envir = .GlobalEnv)
  assign("directory_path", directory_path, envir = .GlobalEnv)
}




create_heatmap <- function(Matrix,
                           results,
                           ALTI,
                           BASSI,
                           n_genes_heat,
                           Metadata_individual,
                           directory2,
                           directory_path,
                           Cell_type,
                           args) {
  
  if (length(c(ALTI$gene, BASSI$gene)) > 2) {
    
    # Select top genes if needed
    if (nrow(ALTI) > n_genes_heat) {
      ALTI <- ALTI[1:n_genes_heat, ]
    }
    if (nrow(BASSI) > n_genes_heat) {
      BASSI <- BASSI[1:n_genes_heat, ]
    }
    
    # Prepare heatmap matrix
    Heat <- Matrix[rownames(Matrix) %in% c(ALTI$CpG_island, BASSI$CpG_island), ]
    Cpg_isl <- rownames(Heat)
    Heat <- apply(Heat, 2, as.numeric)
    rownames(Heat) <- Cpg_isl
    
    SAVED <- NULL
    for (ix in Cpg_isl) {
      iu <- grep(ix, results$CpG_island)
      Gene <- as.character(results$gene[iu])
      SAVED <- append(Gene, SAVED)
    }
    SAVED <- as.character(SAVED)
    SAVED <- rev(SAVED)
    samples <- nrow(Heat)
    
    for (ix in 1:samples) {
      if (SAVED[ix] != "") {
        rownames(Heat)[ix] <- SAVED[ix]
      }
    }
    
    # Plot heatmap
    library(ComplexHeatmap)
    setwd(directory2)
    pdf(file = paste0("Heatmap_top_genes_", Cell_type, "_", args[1], "_vs_", args[2], ".pdf"))
    
    library(circlize)
    col_fun <- colorRamp2(c(min(Heat), mean(colMeans(Heat)), max(Heat)),
                          c("blue", "black", "yellow"))
    col_fun(seq(-3, 3))
    
    t <- c("CTRL" = "blue", "SLE" = "red")
    attr(t, "names")[1] <- args[2]
    attr(t, "names")[2] <- args[1]
    
    ha <- HeatmapAnnotation(condition = Metadata_individual$CONDITION,
                            col = list(condition = t))
    p <- Heatmap(Heat, km = 2, name = "SD_score", col = col_fun,
                 clustering_distance_rows = "pearson",
                 clustering_method_rows = "complete",
                 clustering_method_columns = "ward.D",
                 row_dend_width = unit(0.5, "cm"),
                 column_dend_height = unit(60, "mm"),
                 column_names_gp = grid::gpar(fontsize = 6),
                 row_names_gp = grid::gpar(fontsize = 8),
                 top_annotation = ha)
    print(p)
    
    # Redraw with group names as column names
    colnames(Heat) <- Metadata_individual$CONDITION
    ha <- HeatmapAnnotation(condition = Metadata_individual$CONDITION,
                            col = list(condition = t))
    p <- Heatmap(Heat, km = 2, name = "SD_score", col = col_fun,
                 clustering_distance_rows = "pearson",
                 clustering_method_rows = "complete",
                 clustering_method_columns = "ward.D",
                 row_dend_width = unit(0.5, "cm"),
                 column_dend_height = unit(60, "mm"),
                 column_names_gp = grid::gpar(fontsize = 6),
                 row_names_gp = grid::gpar(fontsize = 8),
                 top_annotation = ha)
    print(p)
    
    dev.off()
    
    # Export variables to global environment for downstream use
    assign("ALTI", ALTI, envir = .GlobalEnv)
    assign("BASSI", BASSI, envir = .GlobalEnv)
    assign("Heat", Heat, envir = .GlobalEnv)
    assign("SAVED", SAVED, envir = .GlobalEnv)
  }
  
  # Create PDF for downstream pathway analysis
  pdf(file = paste0(directory_path, "/Pathways_DMG_EnrichR.pdf"), width = 20, height = 9)
}




perform_pathway_analysis <- function(ALTI,
                                     BASSI,
                                     directory_path,
                                     dbs = c("Reactome_2022",
                                             "GO_Biological_Process_2023",
                                             "CODE_and_ChEA_Consensus_TFs_from_ChIP-X")) {
  library(enrichR)
  websiteLive <- getOption("enrichR.live")
  
  # Create TABLES folder if not existing
  dir.create(file.path(directory_path, "TABLES"), showWarnings = FALSE, recursive = TRUE)
  
  # --- UPREGULATED GENES ---
  if (length(ALTI$gene) != 0) {
    if (websiteLive) {
      enriched <- enrichr(as.character(ALTI$gene), dbs)
    }
    for (paths in 1:length(dbs)) {
      if (nrow(enriched[[paths]]) != 0) {
        if (websiteLive) {
          print(enriched[[paths]])
          write.table(
            x = enriched[[paths]],
            file = file.path(directory_path, "TABLES",
                             paste0(names(enriched)[paths], "_upregulated.tsv")),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
          )
          xx <- plotEnrich(enriched[[paths]], showTerms = 25, numChar = 40,
                           y = "Count", orderBy = "P.value",
                           title = paste0("Enrichment_analysis_upregulated_genes_", dbs[paths]))
          print(xx)
        }
      }
    }
    assign("ALTI_enrichment", enriched, envir = .GlobalEnv)
  }
  
  # --- DOWNREGULATED GENES ---
  if (length(BASSI$gene) != 0) {
    if (websiteLive) {
      enriched <- enrichr(as.character(BASSI$gene), dbs)
    }
    for (paths in 1:length(dbs)) {
      if (nrow(enriched[[paths]]) != 0) {
        if (websiteLive) {
          write.table(
            x = enriched[[paths]],
            file = file.path(directory_path, "TABLES",
                             paste0(names(enriched)[paths], "_downregulated.tsv")),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
          )
          xx <- plotEnrich(enriched[[paths]], showTerms = 25, numChar = 40,
                           y = "Count", orderBy = "P.value",
                           title = paste0("Enrichment_analysis_downregulated_genes_", dbs[paths]))
          print(xx)
        }
      }
    }
    assign("BASSI_enrichment", enriched, envir = .GlobalEnv)
  }
  
  dev.off()
}

