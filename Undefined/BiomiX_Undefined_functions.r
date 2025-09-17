process_metadata_matrix <- function(matrix, Metadata_total) {
  # --- Original code, unchanged in logic ---
  Metadata_Bcell <- Metadata_total
  Metadata_Bcell <- Metadata_Bcell[order(Metadata_Bcell$ID), ]
  
  matrixi <- matrix[, 2:ncol(matrix)]
  matrix[, 2:ncol(matrix)] <- matrixi[, order(colnames(matrixi))]
  matrix <- as.data.frame(matrix)
  
  Metadata_Bcell <- Metadata_total
  num <- which(matrix$ID %in% Metadata_Bcell$ID)
  Identifier <- matrix$ID
  matrix <- as.matrix(matrix[num, ])
  
  print(Metadata_Bcell$ID)
  
  num <- which(Metadata_Bcell$ID %in% matrix[, 1])
  Metadata_Bcell <- Metadata_Bcell[num, ]
  
  # If multiple select only the first column
  Metadata_Bcell <- Metadata_Bcell[!duplicated(Metadata_Bcell$ID), ]
  
  Metadata_Bcell <- Metadata_Bcell[order(Metadata_Bcell$ID), ]
  matrix <- matrix[order(matrix[, 1]), ]
  
  # --- Assign results to global environment ---
  assign("Metadata_Bcell", Metadata_Bcell, envir = .GlobalEnv)
  assign("matrix", matrix, envir = .GlobalEnv)
  assign("Identifier", Identifier, envir = .GlobalEnv)
}



apply_metadata_filter <- function(Metadata, matrix, COMMAND_ADVANCED) {
  # --- Original code, unchanged ---
  METADATA_FILT <- !is.na(COMMAND_ADVANCED[3, grep("*.FILTERING.*", colnames(COMMAND_ADVANCED))])
  METADATA_FILT_INDEX <- grep("*.FILTERING.*", colnames(COMMAND_ADVANCED))
  
  repetition <- 0
  for (meta_filter in METADATA_FILT_INDEX) {
    repetition <- repetition + 1
    
    if (!is.na(COMMAND_ADVANCED[3, grep("*.FILTERING.*", colnames(COMMAND_ADVANCED))])[repetition]) {
      COLNAME <- as.character(COMMAND_ADVANCED[1, meta_filter])
      
      if (as.character(COMMAND_ADVANCED[2, meta_filter]) == "Numerical") {
        To_filter <- as.numeric(unlist(Metadata[, COLNAME]))
        simbol <- substr(as.character(COMMAND_ADVANCED[3, meta_filter]), 1, 1)
        characters_to_remove <- c(">", "<", "=", " ")
        value_threshold <- as.numeric(gsub(paste(characters_to_remove, collapse = "|"), "", 
                                           as.character(COMMAND_ADVANCED[3, meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "<" = function(a, b) a < b,
                                      ">" = function(a, b) a > b,
                                      "=" = function(a, b) a == b,
                                      ">=" = function(a, b) a >= b,
                                      "<=" = function(a, b) a <= b,
                                      NA)
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold), ]
        matrix <- matrix[, comparison_operator(To_filter, value_threshold)]
        
      } else if (as.character(COMMAND_ADVANCED[2, meta_filter]) == "Factors") {
        To_filter <- as.character(unlist(Metadata[, COLNAME]))
        simbol <- substr(as.character(COMMAND_ADVANCED[3, meta_filter]), 1, 2)
        characters_to_remove <- c("!=", "==", " ")
        value_threshold <- as.character(gsub(paste(characters_to_remove, collapse = "|"), "", 
                                             as.character(COMMAND_ADVANCED[3, meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "==" = function(a, b) a == b,
                                      "!=" = function(a, b) a != b,
                                      NA)
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold), ]
        if (!is.null(matrix)) {
          matrix <- matrix[, comparison_operator(To_filter, value_threshold)]
        }
      }
      
    } else {
      print("skip_filtering")
    }
  }
  
  # --- Assign results to global environment ---
  assign("Metadata", Metadata, envir = .GlobalEnv)
  assign("matrix", matrix, envir = .GlobalEnv)
  assign("METADATA_FILT", METADATA_FILT, envir = .GlobalEnv)
  assign("METADATA_FILT_INDEX", METADATA_FILT_INDEX, envir = .GlobalEnv)
}



run_qc_preview <- function(matrix, Metadata, directory) {
  # --- Original code preserved ---
  
  # Add the ID to the first column
  sample_preview <- matrix$ID
  numeric_data <- matrix[, -1]
  
  numeric_data <- apply(numeric_data, 2, as.numeric)
  rownames(numeric_data) <- sample_preview
  
  metadata <- Metadata
  numeric_data <- as.data.frame(numeric_data)
  colnames(metadata)[colnames(metadata) == "CONDITION"] <- "condition"
  
  # Load the QC preview Shiny app
  source(paste(directory, '/BiomiX_preview.r', sep = ""))
  
  browser_analysis <- readLines(paste(directory, '/_INSTALL/CHOISE_BROWSER_pre-view', sep = ""), n = 1)
  
  # Run the Shiny app
  options(browser = get_default_browser())
  print("Pre QC data visualization. If the browser does not open automatically copy and paste the link on a browser")
  print("One completed the analysis modification click on close app to continue the analysis")
  print("ATTENTION!: IF the browser does not open set up the path to your browser on the CHOISE_BROWSER_pre-view file in the BiomiX _INSTALL folder")
  
  Preview <- shiny::runApp(runShinyApp(numeric_data, metadata), launch.browser = TRUE)
  
  # Update matrix and Metadata with the preview results
  matrix <- Preview$matrix
  Metadata <- Preview$metadata
  matrix <- cbind(ID = Metadata$ID, matrix)
  colnames(Metadata)[colnames(Metadata) == "condition"] <- "CONDITION"
  
  # --- Assign results back to global environment ---
  assign("matrix", matrix, envir = .GlobalEnv)
  assign("Metadata", Metadata, envir = .GlobalEnv)
  assign("Preview", Preview, envir = .GlobalEnv)
}







# #### FUNCTION DEFINITION FOR CONTROL AND TESTED SAMPLES----


DMS_REFERENCE <- function(Condition1){
  print(Condition1)
  CON <- matrix[Metadata$CONDITION == Condition1,] #Example CONTROL vs SLE patients
  samp_CON<-row.names(CON) 
  CON <-CON[,c(-1,-2)] 
  CON <-apply(CON,2,as.numeric)
  
  
  
  y=NULL
  for (i in seq(1:length(colnames(CON)))) {
    if (all_same(CON[,i])){
      y <-append(y, 1)
    }else{
      res <- (shapiro.test(CON[,i]))
      y <-append(y, res[["p.value"]])
    }}
  
  CON<-as.data.frame(CON)
  shapiro_test_CON <- y
  shapiro <-shapiro_test_CON < 0.05
  print(paste(length(shapiro[shapiro== TRUE]), "variables out of", length(shapiro), "not having a normal distribution in", args[2]))
  if (length(shapiro[shapiro== TRUE]) > length(shapiro)/2 ) {
    print(paste("Suggested Log normalization before the MOFA analysis for", args[2]))
  }
  
  rownames(CON) <-samp_CON
  return(CON)
}





DMS_TESTED <- function(Condition2){
  print(Condition2)
  SLE <- matrix[Metadata$CONDITION == Condition2,] 
  samp_SLE<-row.names(SLE) 
  SLE <-SLE[,c(-1,-2)]
  SLE <-apply(SLE,2,as.numeric) 
  
  y=NULL
  for (i in seq(1:length(colnames(SLE)))) {
    if (all_same(SLE[,i])){
      y <-append(y, 1)
    }else{
      res <- (shapiro.test(SLE[,i]))
      y <-append(y, res[["p.value"]])
    }}
  
  SLE<-as.data.frame(SLE)
  shapiro_test_SLE <- y
  shapiro <-shapiro_test_SLE < 0.05
  length(shapiro[shapiro== TRUE]) 
  print(paste(length(shapiro[shapiro== TRUE]), "variables out of", length(shapiro), "not having a normal distribution in",args[1]))
  if (length(shapiro[shapiro== TRUE]) > length(shapiro)/2 ) {
    print(paste("Suggested Log normalization before the MOFA analysis for", args[1]))
  }
  #If the p.value is bigger than 0.05 we don't have difference
  #between our and a normal distribution (generated randomly)
  
  rownames(SLE) <-samp_SLE
  
  # Wilcox test CON vs SLE
  pval=NULL
  for (i in 1:ncol(SLE)) {
    res <-wilcox.test(CON[,i],SLE[,i], alternative = "two.sided") 
    pval <-append(pval, res[["p.value"]])
  }
  
  # T test CON vs SLE
  pval_t=NULL
  for (i in 1:ncol(SLE)) {
    res_t <-t.test(CON[,i],SLE[,i], alternative = "two.sided") 
    pval_t <-append(pval_t, res[["p.value"]])
  }
  
  
  
  fold=NULL
  for (i in 1:ncol(CON)) {
    up <-CON[,i]
    up <-up[!is.na(up)]
    down <-SLE[,i]
    down <-down[!is.na(down)]
    FC <- log2(abs(median(down) / median(up)))
    fold <-append(fold, FC)}
  
  
  SLE<-as.data.frame(t(SLE)) #RUN THESE LINES FOR THE CTRL VS SLE COMPARISON
  colnames(SLE) <- samp_SLE
  SLE$shapiro_pvalue <- shapiro_test_SLE
  SLE$p_val_wilcox <- pval
  SLE$p_val_t_test <- pval_t
  SLE$log2FC <- fold
  
  return(SLE)
}


#Function to check if all the variable elements are the same
all_same <- function(x) {
  x<-x[!is.na(x)]
  all(x == x[1])
}






run_differential_analysis <- function(TEST, matrix, Metadata, directory, Cell_type, args, padju, LogFC) {
  # --- Directory and file preparation ---
  dir.create(path = paste(directory, "/Undefined/OUTPUT/", sep = ""), showWarnings = TRUE, recursive = TRUE, mode = "0777")
  dir.create(path = paste(directory, "/Undefined/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep = ""), 
             showWarnings = TRUE, recursive = TRUE, mode = "0777")
  directory2 <- paste(directory, "/Undefined/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep = "")
  
  setwd(directory2)
  
  # --- Statistical calculations ---
  total <- TEST
  total$padj_wilcox <- p.adjust(total$p_val_wilcox, method = "fdr")
  total$padj_t_test <- p.adjust(total$p_val_t_test, method = "fdr")
  totalshOK <- total %>% dplyr::filter(padj_t_test < padju)
  total <- total %>% dplyr::arrange(padj_t_test)
  x <- colnames(total) %in% matrix$ID
  total <- total[, !x]
  total_3 <- total
  total_3 <- total_3[, -2:-3]
  write.table(total_3, paste(directory2, "/", Cell_type, "_", args[1], "_vs_", args[2], "_results.tsv", sep = ""),
              quote = FALSE, row.names = FALSE, sep = "\t")
  
  total$padj_wilcox <- as.numeric(total$padj_wilcox)
  total$log2FC <- as.numeric(total$log2FC)
  
  # --- Identify significantly up/downregulated variables ---
  ALTI <- subset(total, padj_wilcox < padju & log2FC > LogFC)
  ALTI <- ALTI[order(ALTI$padj_wilcox), ]
  
  BASSI <- subset(total, padj_wilcox < padju & log2FC < -LogFC)
  BASSI <- BASSI[order(BASSI$padj_wilcox), ]
  
  x <- total$NAME %in% ALTI$NAME
  NO <- total[!x, ]
  x <- NO$NAME %in% BASSI$NAME
  NO <- NO[!x, ]
  
  # --- Write results tables ---
  write.table(ALTI, file = paste(Cell_type, "_", args[1], "_", args[2], "_VARIABLES_UP.tsv", sep = ""), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(BASSI, file = paste(Cell_type, "_", args[1], "_", args[2], "_VARIABLES_DOWN.tsv", sep = ""), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE), c("ID", "CONDITION")],
              file = paste(args[1], "_", args[2], "_", Cell_type, "_metadata.tsv", sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  SKIP_ALTI <- nrow(ALTI) > 0
  SKIP_BASSI <- nrow(BASSI) > 0
  
  pdf(file = paste("plot_Undefined_", args[1], "_", args[2], "_", Cell_type, ".pdf", sep = ""))
  
  library(ggplot2)
  library(ggrepel)
  
  # --- Volcano Plot for Wilcox Test ---
  p <- ggplot() +
    ggtitle(paste("Variables analysis", Cell_type, "_", args[1], "_vs_", args[2], sep = "")) +
    theme(
      plot.title = element_text(color = "black", size = 14, face = "bold.italic", hjust = 0.5),
      axis.title.x = element_text(color = "black", size = 14, face = "bold"),
      axis.title.y = element_text(color = "black", size = 14, face = "bold")
    ) +
    ylab("-log10(padj_wilcox)") + xlab("log2FC")
  
  if (SKIP_ALTI) {
    p <- p + geom_point(data = ALTI, aes(x = log2FC, y = -log10(padj_wilcox)), color = "red") +
      geom_text_repel(data = ALTI[1:min(25, nrow(ALTI)), ], aes(x = log2FC, y = -log10(padj_wilcox), label = NAME),
                      hjust = 0, vjust = 0, size = 3, force = 10, force_pull = 1, max.overlaps = 5)
  }
  if (SKIP_BASSI) {
    p <- p + geom_point(data = BASSI, aes(x = log2FC, y = -log10(padj_wilcox)), color = "blue") +
      geom_text_repel(data = BASSI[1:min(25, nrow(BASSI)), ], aes(x = log2FC, y = -log10(padj_wilcox), label = NAME),
                      hjust = 0, vjust = 0, size = 3, force = 20, force_pull = 1, max.overlaps = 5)
  }
  
  p <- p + geom_point(data = NO, aes(x = log2FC, y = -log10(padj_wilcox)), color = "black")
  p2 <- p + geom_vline(xintercept = c(-LogFC, LogFC), col = "red") +
    geom_hline(yintercept = -log10(padju), col = "red")
  
  print(p2)
  
  # --- Volcano Plot for T-test ---
  ALTI_2 <- subset(total, padj_t_test < padju & log2FC > LogFC)
  ALTI_2 <- ALTI_2[order(ALTI_2$padj_t_test), ]
  
  BASSI_2 <- subset(total, padj_t_test < padju & log2FC < -LogFC)
  BASSI_2 <- BASSI_2[order(BASSI_2$padj_t_test), ]
  
  x <- total$NAME %in% ALTI_2$NAME
  NO_2 <- total[!x, ]
  x <- NO_2$NAME %in% BASSI_2$NAME
  NO_2 <- NO_2[!x, ]
  
  SKIP_ALTI_2 <- nrow(ALTI_2) > 0
  SKIP_BASSI_2 <- nrow(BASSI_2) > 0
  
  p <- ggplot() +
    ggtitle(paste("Variables analysis", Cell_type, "_", args[1], "_vs_", args[2], sep = "")) +
    theme(
      plot.title = element_text(color = "black", size = 14, face = "bold.italic", hjust = 0.5),
      axis.title.x = element_text(color = "black", size = 14, face = "bold"),
      axis.title.y = element_text(color = "black", size = 14, face = "bold")
    ) +
    ylab("-log10(padj_t_test)") + xlab("log2FC")
  
  if (SKIP_ALTI_2) {
    p <- p + geom_point(data = ALTI_2, aes(x = log2FC, y = -log10(padj_t_test)), color = "red") +
      geom_text_repel(data = ALTI_2[1:min(25, nrow(ALTI_2)), ], aes(x = log2FC, y = -log10(padj_t_test), label = NAME),
                      hjust = 0, vjust = 0, size = 3, force = 10, force_pull = 1, max.overlaps = 5)
  }
  if (SKIP_BASSI_2) {
    p <- p + geom_point(data = BASSI_2, aes(x = log2FC, y = -log10(padj_t_test)), color = "blue") +
      geom_text_repel(data = BASSI_2[1:min(25, nrow(BASSI_2)), ], aes(x = log2FC, y = -log10(padj_t_test), label = NAME),
                      hjust = 0, vjust = 0, size = 3, force = 20, force_pull = 1, max.overlaps = 5)
  }
  
  p <- p + geom_point(data = NO_2, aes(x = log2FC, y = -log10(padj_t_test)), color = "black")
  p2 <- p + geom_vline(xintercept = c(-LogFC, LogFC), col = "red") +
    geom_hline(yintercept = -log10(padju), col = "red")
  
  print(p2)
  
  dev.off()
  
  # --- Assign variables to global environment ---
  assign("total", total, envir = .GlobalEnv)
  assign("ALTI", ALTI, envir = .GlobalEnv)
  assign("BASSI", BASSI, envir = .GlobalEnv)
  assign("NO", NO, envir = .GlobalEnv)
  assign("ALTI_2", ALTI_2, envir = .GlobalEnv)
  assign("BASSI_2", BASSI_2, envir = .GlobalEnv)
  assign("NO_2", NO_2, envir = .GlobalEnv)
  assign("directory2", directory2, envir = .GlobalEnv)
}



generate_heatmap_variables <- function(ALTI, BASSI, matrix, Metadata, directory2, Cell_type, args, Heatmap_genes) {
  ### HEATMAP ###
  if (length(c(ALTI$NAME, BASSI$NAME)) > 2) {
    
    # Limit number of genes included in heatmap
    if (nrow(ALTI) > Heatmap_genes) {
      ALTI <- ALTI[1:Heatmap_genes, ]
    }
    if (nrow(BASSI) > Heatmap_genes) {
      BASSI <- BASSI[1:Heatmap_genes, ]
    }
    
    # Heatmap input: select only columns corresponding to ALTI/BASSI genes
    Heat <- matrix[, colnames(matrix) %in% c(ALTI$NAME, BASSI$NAME)]
    samples <- rownames(Heat)
    Heat <- apply(Heat, 2, as.numeric)
    rownames(Heat) <- samples
    Heat <- t(Heat)
    Heat[Heat == 0] <- 1
    Heat <- apply(Heat, 2, log10)
    
    library(ComplexHeatmap)
    library(circlize)
    
    setwd(directory2)
    pdf(file = paste("Heatmap_top_variables_", Cell_type, "_", args[2], "_vs_", args[1], ".pdf", sep = ""))
    
    # Color scale function
    col_fun <- colorRamp2(
      c(min(Heat, na.rm = TRUE), mean(colMeans(Heat, na.rm = TRUE)), max(Heat, na.rm = TRUE)),
      c("blue", "black", "yellow")
    )
    col_fun(seq(-3, 3))
    
    # Conditions color annotation
    t <- c("CTRL" = "blue", "SLE" = "red")
    attr(t, "names")[1] <- args[2]
    attr(t, "names")[2] <- args[1]
    
    ha <- HeatmapAnnotation(condition = Metadata$CONDITION,
                            col = list(condition = t))
    
    Heat[is.na(Heat)] <- 0 # Replace NA with 0 to avoid plotting errors
    
    # Generate heatmap
    p <- Heatmap(
      Heat, km = 2, name = "SD_score", col = col_fun,
      clustering_distance_rows = "pearson",
      clustering_method_rows = "complete",
      clustering_method_columns = "ward.D",
      row_dend_width = unit(0.5, "cm"),
      column_dend_height = unit(60, "mm"),
      column_names_gp = grid::gpar(fontsize = 6),
      row_names_gp = grid::gpar(fontsize = 8),
      top_annotation = ha
    )
    print(p)
    
    # Add group visualization
    colnames(Heat) <- Metadata$CONDITION
    ha <- HeatmapAnnotation(condition = Metadata$CONDITION,
                            col = list(condition = t))
    p <- Heatmap(
      Heat, km = 2, name = "SD_score", col = col_fun,
      clustering_distance_rows = "pearson",
      clustering_method_rows = "complete",
      clustering_method_columns = "ward.D",
      row_dend_width = unit(0.5, "cm"),
      column_dend_height = unit(60, "mm"),
      column_names_gp = grid::gpar(fontsize = 6),
      row_names_gp = grid::gpar(fontsize = 8),
      top_annotation = ha
    )
    print(p)
    
    dev.off()
    
    # Assign objects back to global environment (so updates propagate)
    assign("ALTI", ALTI, envir = .GlobalEnv)
    assign("BASSI", BASSI, envir = .GlobalEnv)
    assign("Heat", Heat, envir = .GlobalEnv)
  }
}

