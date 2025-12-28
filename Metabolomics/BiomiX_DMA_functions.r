#Functions metabolomics pipeline BiomiX_DMA.r

#LOAD ANNOTATION FUNCTION

load_annotation <- function(ANNOTATION, i, COMMAND_ADVANCED) {
  # Function to load MS1 or MS2 annotation files and prepare `annotation`
  
  if (ANNOTATION == "MS1") {
    input_ms1 <- which(substr(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX, 1, 1) %in% i)
    if (length(input_ms1) != 0) {
      if (file.exists(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES[input_ms1])) {
        
        if (grepl("\\.xlsx$|\\.xls$", COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES[input_ms1])) {
          annotation <- readxl::read_excel(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES[input_ms1])
          print("Annotation ms1 Excel File read successfully!")
        } else {
          annotation <- vroom::vroom(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES[input_ms1], delim = "\t", col_names = TRUE)
        }
        
        if (ncol(annotation) == 3) {
          if (colnames(annotation)[3] == "RT_sec") {
            annotation$RT_min <- as.numeric(annotation$RT_sec) / 60
          }
          if (colnames(annotation)[3] == "RT_min") {
            annotation$RT_sec <- as.numeric(annotation$RT_min) * 60
          }
        }
      }
    }
  }
  
  if (ANNOTATION == "MS2") {
    input_ms1 <- which(substr(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX, 1, 1) %in% i)
    if (length(input_ms1) != 0) {
      if (file.exists(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2[input_ms1])) {
        
        if (grepl("\\.xlsx$|\\.xls$", COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2[input_ms1])) {
          annotation <- readxl::read_excel(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2[input_ms1])
          print("Annotation ms1 Excel File read successfully!")
        } else {
          annotation <- vroom::vroom(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2[input_ms1], delim = "\t", col_names = TRUE)
        }
        
        if (ncol(annotation) == 3) {
          if (colnames(annotation)[3] == "RT_sec") {
            annotation$RT_min <- as.numeric(annotation$RT_sec) / 60
          }
          if (colnames(annotation)[3] == "RT_min") {
            annotation$RT_sec <- as.numeric(annotation$RT_min) * 60
          }
        }
      }
    }
  }
  
  # Assign results to global environment to preserve behavior
  assign("annotation", annotation, envir = .GlobalEnv)
  assign("input_ms1", input_ms1, envir = .GlobalEnv)
}









#HMDB METABOLITE FILTERING 

load_biofluid_metabolites <- function(label, biofluid_name, local_filename, url, variable_name, directory2) {
  # Checks if a biofluid is mentioned in the label, then loads (or downloads) the metabolite file.
  
  if (stringr::str_detect(label, stringr::fixed(biofluid_name, ignore_case = TRUE))) {
    if (file.exists(local_filename)) {
      print("File available locally, using the local version")
      metabolite_data <- vroom::vroom(file.path(directory2, local_filename), delim = "\t", col_names = TRUE)
    } else {
      print("File unavailable locally, downloading it from HMDB database")
      options(timeout = 6000)
      metabolite_data <- vroom::vroom(url(url), delim = ",", col_names = TRUE)
      metabolite_data <- as.data.frame(metabolite_data)
      write.table(metabolite_data, file.path(directory2, local_filename), quote = FALSE, row.names = FALSE, sep = "\t")
    }
    # Assign to global environment with dynamic name
    assign(variable_name, metabolite_data, envir = .GlobalEnv)
  }
}




#Function for sample Selection by label

process_sample_selection <- function(selection_samples, Cell_type, Metadata_total, matrix) {
  if (selection_samples == "YES") {
    num <- grep(Cell_type, Metadata_total$ID_CELL_TYPE)
    Metadata_Bcell <- Metadata_total[num, ]
    Metadata_total <- Metadata_total[num, ]
    Identifier <- matrix$ID
    num <- which(colnames(matrix) %in% Metadata_Bcell$ID_CELL_TYPE)
    matrix <- as.matrix(matrix[, num])
    print(Metadata_Bcell$ID_CELL_TYPE)
    num <- which(Metadata_Bcell$ID_CELL_TYPE %in% colnames(matrix))
    Metadata_Bcell <- Metadata_Bcell[num, ]
    
    # If multiple select only the first column
    Metadata_Bcell <- Metadata_Bcell[!duplicated(Metadata_Bcell$ID), ]
    
    Metadata_Bcell <- Metadata_Bcell[order(Metadata_Bcell$ID), ]
    matrix <- matrix[, order(colnames(matrix))]
    
  } else {
    print("No samples selection")
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
  }
  
  # Assign results to global environment
  assign("Metadata_total", Metadata_total, envir = .GlobalEnv)
  assign("Metadata_Bcell", Metadata_Bcell, envir = .GlobalEnv)
  assign("matrix", matrix, envir = .GlobalEnv)
  assign("Identifier", Identifier, envir = .GlobalEnv)
}




# FILTERING SAMPLES BASED ON METADATA CRITERIA SELECTED
# Defined in advanced options (metadata section)

filter_samples_by_metadata <- function(COMMAND_ADVANCED, Metadata, matrix) {

  
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
                                      NA
        )
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold), ]
        matrix <- matrix[comparison_operator(To_filter, value_threshold), ]
        
      } else if (as.character(COMMAND_ADVANCED[2, meta_filter]) == "Factors") {
        To_filter <- as.character(unlist(Metadata[, COLNAME]))
        simbol <- substr(as.character(COMMAND_ADVANCED[3, meta_filter]), 1, 2)
        characters_to_remove <- c("!=", "==", " ")
        value_threshold <- as.character(gsub(paste(characters_to_remove, collapse = "|"), "", 
                                             as.character(COMMAND_ADVANCED[3, meta_filter])))
        
        comparison_operator <- switch(simbol,
                                      "==" = function(a, b) a == b,
                                      "!=" = function(a, b) a != b,
                                      NA
        )
        
        Metadata <- Metadata[comparison_operator(To_filter, value_threshold), ]
        matrix <- matrix[comparison_operator(To_filter, value_threshold), ]
      }
      
    } else {
      print("skip_filtering")
    }
  }
  
  # Assign updated objects back to the global environment
  assign("Metadata", Metadata, envir = .GlobalEnv)
  assign("matrix", matrix, envir = .GlobalEnv)
}





#FUNCTION QC VISUALISATION

run_qc_preview <- function(COMMAND, i, matrix, Metadata, directory) {
  if (COMMAND$PREVIEW[i] == "YES") {
    
    # Block to open the Quality control windows to select filtering, normalization and visualize the QC.
    # ADD the ID to the first column
    sample_preview <- matrix$ID
    numeric_data <- matrix[, -1]
    
    numeric_data <- apply(numeric_data, 2, as.numeric)
    rownames(numeric_data) <- sample_preview
    
    metadata <- Metadata
    numeric_data <- as.data.frame(numeric_data)
    colnames(metadata)[colnames(metadata) == "CONDITION"] <- "condition"
    
    # Source the BiomiX_preview script to load the runShinyApp() function
    source(paste(directory, '/BiomiX_preview.r', sep = ""))
    
    browser_analysis <- readLines(paste(directory, '/_INSTALL/CHOISE_BROWSER_pre-view', sep = ""), n = 1)
    
    # Call the runShinyApp function with numeric_data and metadata
    options(browser = get_default_browser())
    print("Pre QC data visualization. If the browser does not open automatically copy and paste the link on a browser")
    print("Once completed the analysis modification click on close app to continue the analysis")
    print("ATTENTION!: IF the browser does not open set up the path to your browser on the CHOISE_BROWSER_pre-view file in the BiomiX _INSTALL folder")
    
    Preview <- shiny::runApp(runShinyApp(numeric_data, metadata), launch.browser = TRUE)
    
    matrix <- Preview$matrix
    Metadata <- Preview$metadata
    
    Metadata<-Metadata[!Metadata$condition == "QC",]
    matrix <-matrix[!Metadata$condition == "QC",]
    
    matrix <- cbind(ID = Metadata$ID, matrix)
    colnames(Metadata)[colnames(Metadata) == "condition"] <- "CONDITION"
    
  } else {
    print("no QC pre-visualization")
  }
  
  # Assign back to global environment
  assign("matrix", matrix, envir = .GlobalEnv)
  assign("Metadata", Metadata, envir = .GlobalEnv)
}




DMS_REFERENCE <- function(Condition1){
  print(Condition1)
  CON <- matrix[Metadata$CONDITION == Condition1,] #Example CONTROL vs SLE patients
  samp_CON<-row.names(CON) #RUN IT ALWAYS + ONE DISEASE
  CON <-CON[,c(-1,-2)] #RUN IT ALWAYS + ONE DISEASE
  CON <-apply(CON,2,as.numeric) #RUN IT ALWAYS + ONE DISEASE
  
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
  SLE <- matrix[Metadata$CONDITION == Condition2,] #RUN LINE 58 and 59 for CTRL vs SLE
  samp_SLE<-row.names(SLE) #RUN LINE 80 and 81 for CTRL vs SLE
  SLE <-SLE[,c(-1,-2)] #RUN LINE 109 and 110 for CTRL vs SLE
  SLE <-apply(SLE,2,as.numeric) #RUN LINE 130 and 131 for CTRL vs SLE
  
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
  SLE$p_val <- pval
  SLE$log2FC <- fold
  
  return(SLE)
}


#Function to check if all the variable elements are the same
all_same <- function(x) {
  x<-x[!is.na(x)]
  all(x == x[1])
}






# DGE FUNCTION FOR ANNOTATED PIPELINE


run_dge_analysis_annotated <- function(TEST, matrix, Metadata, Cell_type, args, directory, padju, LogFC) {
  # Create output directories
  dir.create(path = paste0(directory, "/Metabolomics/OUTPUT/"), showWarnings = TRUE, recursive = TRUE, mode = "0777")
  dir.create(path = paste0(directory, "/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2]),
             showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  directory2 <- paste0(directory, "/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2])
  setwd(directory2)
  
  total <- TEST
  total$padj <- p.adjust(total$p_val, method = "fdr")
  totalshOK <- dplyr::filter(total, padj < padju)
  total <- dplyr::arrange(total, padj)
  x <- colnames(total) %in% matrix$ID
  total <- total[, !x]
  
  
  
  #Inclusion HMDB
  total <- total %>%
    left_join(Mart_metabolome ,
              by = c("NAME" = "HMDB")) %>%
    mutate(name = if_else(!is.na(name), name, NAME))
  
  colnames(total)[which(colnames(total) == "name")] <- "Name"

  
  total_3 <- total
  colnames(total_3)[colnames(total_3) == "NAME"] <- "NAME.x"
  
  write.table(total_3, paste0(directory2, "/", Cell_type, "_", args[1], "_vs_", args[2], "_results.tsv"),
              quote = FALSE, row.names = FALSE, sep = "\t")
  
  total$padj <- as.numeric(total$padj)
  total$log2FC <- as.numeric(total$log2FC)
  
  
  total_min <- total[, c("NAME", "log2FC", "p_val", "padj")]
  total_min <- dplyr::arrange(total_min, padj)
  write.table(total_min, paste0(directory2, "/", Cell_type, "_", args[1], "_vs_", args[2], "_peak_statistics.tsv"),
              quote = FALSE, row.names = FALSE, sep = "\t")
  
  # Identify significant hits
  ALTI <- subset(total, padj < padju & log2FC > LogFC)
  ALTI <- ALTI[order(ALTI$padj), ]
  
  BASSI <- subset(total, padj < padju & log2FC < -LogFC)
  BASSI <- BASSI[order(BASSI$padj), ]
  
  NO <- total[!(total$NAME %in% c(ALTI$NAME, BASSI$NAME)), ]
  
  # Write results
  write.table(ALTI$Name, paste0(Cell_type, "_", args[1], "_", args[2], "_UP.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(BASSI$Name, paste0(Cell_type, "_", args[1], "_", args[2], "_DOWN.tsv"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(ALTI, paste0(Cell_type, "_", args[1], "_", args[2], "_METABUP.tsv"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(BASSI, paste0(Cell_type, "_", args[1], "_", args[2], "_METABDOWN.tsv"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE), c("ID", "CONDITION")],
              file = paste0(args[1], "_", args[2], "_", Cell_type, "_metadata.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Flags
  SKIP_ALTI <- nrow(ALTI) > 0
  SKIP_BASSI <- nrow(BASSI) > 0
  
  # Volcano Plot
  pdf(file = paste0("plot_DMA_", args[1], "_", args[2], "_", Cell_type, ".pdf"))
  library(ggplot2)
  library(ggrepel)
  
  p <- ggplot() +
    ggtitle(paste("DGE", Cell_type, "_", args[1], "_vs_", args[2], sep = "")) +
    theme(plot.title = element_text(color = "black", size = 14, face = "bold.italic", hjust = 0.5),
          axis.title.x = element_text(color = "black", size = 14, face = "bold"),
          axis.title.y = element_text(color = "black", size = 14, face = "bold")) +
    ylab("-log10(padj)") + xlab("log2FC")
  
  if (SKIP_ALTI) {
    p <- p + geom_point(data = ALTI, aes(x = log2FC, y = -log10(padj)), color = "red") +
      geom_text_repel(data = ALTI[1:min(25, nrow(ALTI)), ], aes(x = log2FC, y = -log10(padj), label = Name),
                      hjust = 0, vjust = 0, size = 3, max.overlaps = 5)
  }
  if (SKIP_BASSI) {
    p <- p + geom_point(data = BASSI, aes(x = log2FC, y = -log10(padj)), color = "blue") +
      geom_text_repel(data = BASSI[1:min(25, nrow(BASSI)), ], aes(x = log2FC, y = -log10(padj), label = Name),
                      hjust = 0, vjust = 0, size = 3, max.overlaps = 5)
  }
  
  p <- p + geom_point(data = NO, aes(x = log2FC, y = -log10(padj)), color = "black") +
    geom_vline(xintercept = c(-LogFC, LogFC), col = "red") +
    geom_hline(yintercept = -log10(padju), col = "red")
  
  print(p)
  dev.off()
  
  # Assign results back to global environment
  assign("ALTI", ALTI, envir = .GlobalEnv)
  assign("BASSI", BASSI, envir = .GlobalEnv)
  assign("NO", NO, envir = .GlobalEnv)
  assign("total", total, envir = .GlobalEnv)
  assign("directory2", directory2, envir = .GlobalEnv)
}




#HEATMAP FOR ANNOTATED PIPELINE


generate_heatmap_annotated <- function(ALTI, BASSI, matrix, Metadata, Heatmap_genes, Cell_type, args, directory2) {
  if (length(c(ALTI$NAME, BASSI$NAME)) > 2) {
    # Limit number of top genes
    if (nrow(ALTI) > Heatmap_genes) ALTI <- ALTI[1:Heatmap_genes, ]
    if (nrow(BASSI) > Heatmap_genes) BASSI <- BASSI[1:Heatmap_genes, ]
    
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
    pdf(file = paste0("Heatmap_top_genes_", Cell_type, "_", args[2], "_vs_", args[1], ".pdf"))
    
    col_fun <- colorRamp2(c(min(Heat, na.rm = TRUE),
                            mean(colMeans(Heat, na.rm = TRUE)),
                            max(Heat, na.rm = TRUE)),
                          c("blue", "black", "yellow"))
    
    t <- c("CTRL" = "blue", "SLE" = "red")
    attr(t, "names")[1] <- args[2]
    attr(t, "names")[2] <- args[1]
    
    ha <- HeatmapAnnotation(condition = Metadata$CONDITION, col = list(condition = t))
    Heat[is.na(Heat)] <- 0
    
    p <- Heatmap(Heat, km = 2, name = "SD_score", col = col_fun,
                 clustering_distance_rows = "pearson", clustering_method_rows = "complete",
                 clustering_method_columns = "ward.D",
                 row_dend_width = unit(0.5, "cm"), column_dend_height = unit(60, "mm"),
                 column_names_gp = grid::gpar(fontsize = 6),
                 row_names_gp = grid::gpar(fontsize = 8),
                 top_annotation = ha)
    print(p)
    dev.off()
  }
}



#METHPATH PIPELINE FUNCTION

run_metpath_pipeline_annotated <- function(total, COMMAND_ADVANCED, Cell_type, args, directory2, hmdb_pathway, kegg_hsa_pathway) {
  library(metpath)
  library(tidyverse)
  library(dplyr)
  
  
  
  if (nrow(query_id) != 0) {
    if (COMMAND_ADVANCED[2, 9] == "HMDB") {
      colnames(query_id)[1] <- "HMDB"
      
      dir.create(path = paste0(directory2, "/Pathway_analysis/", "HMDB_", Cell_type, "_", args[1], "_vs_", args[2]),
                 showWarnings = TRUE, recursive = TRUE, mode = "0777")
      setwd(paste0(directory2, "/Pathway_analysis/", "HMDB_", Cell_type, "_", args[1], "_vs_", args[2]))
      
      pathway_class_HMDB <- metpath::pathway_class(hmdb_pathway)
      
      gc()
      remain_idx <- which(unlist(pathway_class_HMDB) == "Metabolic;primary_pathway")
      hmdb_pathway <- hmdb_pathway[remain_idx]
      
      result <- enrich_hmdb(
        query_id = unique(query_id$HMDB),
        query_type = "compound",
        id_type = "HMDB",
        pathway_database = hmdb_pathway,
        only_primary_pathway = TRUE,
        p_cutoff = 0.05,
        p_adjust_method = "BH",
        threads = as.numeric(COMMAND_ADVANCED[3, 3])
      )
      
      if (length(result) != 0) {
        pdf(file = paste0("Pathway_analysis_HMDB", args[1], "_", args[2], "_", Cell_type, ".pdf"))
        
        x <- enrich_bar_plot(object = result, x_axis = "p_value_adjust", cutoff = 1.1, top = 10)
        print(x)
        
        x <- enrich_scatter_plot(object = result)
        print(x)
        
        write.table(result@result,
                    paste0(directory2, "/Pathway_analysis/", "HMDB_", Cell_type, "_", args[1], "_vs_", args[2], "/HMDB_table_results"),
                    quote = FALSE, row.names = FALSE, sep = "\t")
        dev.off()
        gc()
      }
    }
    
    if (COMMAND_ADVANCED[2, 9] == "KEGG") {
      colnames(query_id)[1] <- "KEGG"
      pathway_class_KEGG <- metpath::pathway_class(kegg_hsa_pathway)
      
      dir.create(path = paste0(directory2, "/Pathway_analysis/", "KEGG_", Cell_type, "_", args[1], "_vs_", args[2]),
                 showWarnings = TRUE, recursive = TRUE, mode = "0777")
      setwd(paste0(directory2, "/Pathway_analysis/", "KEGG_", Cell_type, "_", args[1], "_vs_", args[2]))
      
      remain_idx <- pathway_class_KEGG %>%
        unlist() %>%
        stringr::str_detect("Disease") %>%
        `!`() %>%
        which()
      
      pathway_database <- kegg_hsa_pathway[remain_idx]
      
      result <- enrich_kegg(
        query_id = unique(query_id$KEGG),
        query_type = "compound",
        id_type = "KEGG",
        pathway_database = pathway_database,
        p_cutoff = 0.05,
        p_adjust_method = "BH",
        threads = as.numeric(COMMAND_ADVANCED[3, 3])
      )
      
      if (length(result) != 0) {
        pdf(file = paste0("Pathway_analysis_KEGG", args[1], "_", args[2], "_", Cell_type, ".pdf"))
        
        x <- enrich_bar_plot(object = result, x_axis = "p_value_adjust", cutoff = 1.1, top = 10)
        print(x)
        
        x <- enrich_scatter_plot(object = result)
        print(x)
        
        write.table(result@result,
                    paste0(directory2, "/Pathway_analysis/", "KEGG_", Cell_type, "_", args[1], "_vs_", args[2], "/KEGG_table_results.tsv"),
                    quote = FALSE, row.names = FALSE, sep = "\t")
        dev.off()
      }
    }
  } else {
    print("NO STATISTICAL SIGNIFICANT PEAK FOR METABOLITE PATHWAY ANALYSIS")
  }
  
  # Make query_id available globally for the next steps
  assign("query_id", query_id, envir = .GlobalEnv)
}









#METABOANALYST PIPELINE FUNCTION

run_metaboanalyst_pipeline_annotated <- function(query_id, COMMAND_ADVANCED, Cell_type, args, directory, directory2) {
  
#Enrichment analysis
dir.create(path = paste(directory2,"/MetaboAnalyst/", "Enrichment_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
setwd(paste(directory2,"/MetaboAnalyst/", "Enrichment_Analysis", sep =""))

if (COMMAND_ADVANCED[2,9] == "HMDB" ){
  write.table(x= query_id$HMDB[!is.na(query_id$HMDB)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
if (COMMAND_ADVANCED[2,9] == "KEGG" ){
  write.table(x= query_id$Kegg[!is.na(query_id$Kegg )] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
if (COMMAND_ADVANCED[2,9] == "compound_name" ){
  write.table(x= query_id$Name[!is.na(query_id$Name )] , file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}


#Pathway analysis
dir.create(path = paste(directory2,"/MetaboAnalyst/", "Pathway_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
setwd(paste(directory2,"/MetaboAnalyst/", "Pathway_Analysis", sep =""))

if (COMMAND_ADVANCED[2,9] == "HMDB" ){
  write.table(x= query_id$HMDB[!is.na(query_id$HMDB)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
if (COMMAND_ADVANCED[2,9] == "KEGG" ){
  write.table(x= query_id$Kegg[!is.na(query_id$Kegg )] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
if (COMMAND_ADVANCED[2,9] == "compound_name" ){
  write.table(x= query_id$Name[!is.na(query_id$Name )] , file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)}


#Joint-Pathway analysis
dir.create(path = paste(directory2,"/MetaboAnalyst/", "Joint_Pathway_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
setwd(paste(directory2,"/MetaboAnalyst/", "Joint_Pathway_Analysis", sep =""))

query_id_select <- query_id[,c(1,2)]

if (COMMAND_ADVANCED[2,9] == "HMDB" ){
  write.table(x= query_id_select, file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
if (COMMAND_ADVANCED[2,9] == "KEGG" ){
  write.table(x= query_id_select , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
if (COMMAND_ADVANCED[2,9] == "compound_name" ){
  write.table(x= query_id_select, file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}


#LATE INTEGRATION TRANSCRIPTOMICS

regeX <- paste("*",args[1],"_vs_",args[2], sep="")
regeX2 <- paste(args[1],"_vs_",args[2], sep="")
files <- grep(regeX,list.files(paste(directory,"/Transcriptomics/OUTPUT",sep="")),value=TRUE)

if (length(files) != 0){
  for (fil in 1:length(files)){
    print(fil)
    files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2, sep="")),value=TRUE)
    LIST_GENE<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[1], sep=""), delim="\t")
    LIST_GENE2<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[2], sep=""), delim="\t")
    LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
    
    Sources<-str_split(files2[1], "_")[[1]][1]
    print(Sources)
    
    dir.create(path = paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Transcriptomes_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Transcriptomes_availables/", Sources, sep =""))
    
    write.table(x= LIST_GENE3[c("Gene.name", "log2FoldChange")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    
  }}

#LATE INTEGRATION METHYLOMICS

regeX <- paste("*",args[1],"_vs_",args[2], sep="")
regeX2 <- paste(args[1],"_vs_",args[2], sep="")
files <- grep(regeX,list.files(paste(directory,"/Methylomics/OUTPUT",sep="")),value=TRUE)

if (length(files) != 0){
  for (fil in 1:length(files)){
    print(fil)
    files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Methylomics/OUTPUT/", files[fil], sep="")),value=TRUE)
    LIST_GENE<-vroom(paste(directory,"/Methylomics/OUTPUT/", files[fil],"/",files2[1], sep=""), delim="\t")
    LIST_GENE2<-vroom(paste(directory,"/Methylomics/OUTPUT/", files[fil],"/",files2[2], sep=""), delim="\t")
    LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
    
    Sources<-str_split(files2[1], "_")[[1]][1]
    print(Sources)
    
    dir.create(path = paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Methylomics_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Methylomics_availables/", Sources, sep =""))
    
    write.table(x= LIST_GENE3[c("gene", "logFC")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    
  }}


#Network analysis
dir.create(path = paste(directory2,"/MetaboAnalyst/", "Network_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
setwd(paste(directory2,"/MetaboAnalyst/", "Network_Analysis", sep =""))

query_id_select <- query_id[,c(1,2)]

if (COMMAND_ADVANCED[2,9] == "HMDB" ){
  write.table(x= query_id_select, file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}
if (COMMAND_ADVANCED[2,9] == "KEGG" ){
  write.table(x= query_id_select , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
} #Sometimes does work, is not able to write the file, exclusively in the network section (?)
if (COMMAND_ADVANCED[2,9] == "compound_name" ){
  write.table(x= query_id_select, file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)}


#LATE INTEGRATION TRANSCRIPTOMICS

regeX <- paste("*",args[1],"_vs_",args[2], sep="")
regeX2 <- paste(args[1],"_vs_",args[2], sep="")
files <- grep(regeX,list.files(paste(directory,"/Transcriptomics/OUTPUT",sep="")),value=TRUE)

if (length(files) != 0){
  for (fil in 1:length(files)){
    print(fil)
    files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2, sep="")),value=TRUE)
    LIST_GENE<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[1], sep=""), delim="\t")
    LIST_GENE2<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[2], sep=""), delim="\t")
    LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
    
    Sources<-str_split(files2[1], "_")[[1]][1]
    print(Sources)
    
    dir.create(path = paste(directory2,"/MetaboAnalyst/Network_Analysis/Transcriptomes_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/Network_Analysis/Transcriptomes_availables/", Sources, sep =""))
    
    write.table(x= LIST_GENE3[c("Gene.name", "log2FoldChange")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    
  }
}


#LATE INTEGRATION METHYLOMICS

regeX <- paste("*",args[1],"_vs_",args[2], sep="")
regeX2 <- paste(args[1],"_vs_",args[2], sep="")
files <- grep(regeX,list.files(paste(directory,"/Methylomics/OUTPUT",sep="")),value=TRUE)

if (length(files) != 0){
  for (fil in 1:length(files)){
    print(fil)
    files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Methylomics/OUTPUT/", files[fil], sep="")),value=TRUE)
    LIST_GENE<-vroom(paste(directory,"/Methylomics/OUTPUT/", files[fil],"/",files2[1], sep=""), delim="\t")
    LIST_GENE2<-vroom(paste(directory,"/Methylomics/OUTPUT/", files[fil],"/",files2[2], sep=""), delim="\t")
    LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
    
    Sources<-str_split(files2[1], "_")[[1]][1]
    print(Sources)
    
    dir.create(path = paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Methylomics_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Methylomics_availables/", Sources, sep =""))
    
    write.table(x= LIST_GENE3[c("gene", "logFC")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    
  }}

assign("query_id", query_id, envir = .GlobalEnv)
}






# Function for MS/MS Filtering (MS2 PIPELINE)
run_msms_filtering <- function(COMMAND_ADVANCED, directory, MS2_databases) {
  library(metid)
  
  directory3 <- paste(directory, "/Metabolomics", sep = "")
  setwd(directory3)
  
  # Extract parameters
  path <- as.character(COMMAND_ADVANCED[1, 19])
  rt_threshold <- as.numeric(COMMAND_ADVANCED[2, 8])
  
  # Initialize param list
  param <- list()
  
  # Load MS2 databases and create parameter sets
  if (sum(MS2_databases == "HMDB") == 1) {
    HMDA <- load(paste(directory, "/Integration/x_BiomiX_DATABASE/hmdb_database0.0.3.rda", sep = ""))
    param1 <- identify_metabolites_params(
      ms1.ms2.match.mz.tol = as.numeric(COMMAND_ADVANCED[1, 8]),
      # rt.match.tol = rt_threshold,
      ms1.ms2.match.rt.tol = rt_threshold,
      polarity = as.character(COMMAND_ADVANCED[1, 24]),
      ce = "all",
      column = as.character(COMMAND_ADVANCED[3, 23]),
      total.score.tol = 0.5,
      candidate.num = 3,
      threads = as.numeric(COMMAND_ADVANCED[3, 3]),
      database = hmdb_database0.0.3
    )
    param <- append(param, param1)
  }
  
  if (sum(MS2_databases == "MASSBANK") == 1) {
    MASSBANK <- load(paste(directory, "/Integration/x_BiomiX_DATABASE/massbank_database0.0.3.rda", sep = ""))
    param1 <- identify_metabolites_params(
      ms1.ms2.match.mz.tol = as.numeric(COMMAND_ADVANCED[1, 8]),
      # rt.match.tol = rt_threshold,
      ms1.ms2.match.rt.tol = rt_threshold,
      polarity = as.character(COMMAND_ADVANCED[1, 24]),
      ce = "all",
      column = as.character(COMMAND_ADVANCED[3, 23]),
      total.score.tol = 0.5,
      candidate.num = 3,
      threads = as.numeric(COMMAND_ADVANCED[3, 3]),
      database = massbank_database0.0.3
    )
    param <- append(param, param1)
  }
  
  if (sum(MS2_databases == "MONA") == 1) {
    MONA <- load(paste(directory, "/Integration/x_BiomiX_DATABASE/mona_database0.0.3.rda", sep = ""))
    param1 <- identify_metabolites_params(
      ms1.ms2.match.mz.tol = as.numeric(COMMAND_ADVANCED[1, 8]),
      # rt.match.tol = rt_threshold,
      ms1.ms2.match.rt.tol = rt_threshold,
      polarity = as.character(COMMAND_ADVANCED[1, 24]),
      ce = "all",
      column = as.character(COMMAND_ADVANCED[3, 23]),
      total.score.tol = 0.5,
      candidate.num = 3,
      threads = as.numeric(COMMAND_ADVANCED[3, 3]),
      database = mona_database0.0.3
    )
    param <- append(param, param1)
  }
  
  # Make variables available globally
  assign("directory3", directory3, envir = .GlobalEnv)
  assign("path", path, envir = .GlobalEnv)
  assign("rt_threshold", rt_threshold, envir = .GlobalEnv)
  assign("param", param, envir = .GlobalEnv)
}








# ---- Function to Select Peaks and Build Annotation ----
build_msms_annotation <- function(directory, Cell_type, args, param, annotate_result5) {
  
  # Create necessary output directories
  dir.create(path = paste(directory, "/Metabolomics/OUTPUT/", sep = ""), 
             showWarnings = TRUE, recursive = TRUE, mode = "0777")
  dir.create(path = paste(directory, "/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep = ""), 
             showWarnings = TRUE, recursive = TRUE, mode = "0777")
  
  directory2 <- paste(directory, "/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep = "")
  setwd(directory2)
  
  gc()
  unione <- list()
  unione_all <- list()
  iter <- 0
  
  for (dat in 1:length(param)) {
    iter <- iter + 1
    xx <- names(annotate_result5[[dat]]@identification.result)
    saved <- annotate_result5[[dat]]@match.result[
      annotate_result5[[dat]]@match.result$MS2.spectra.name %in% xx, ]
    
    y <- NULL
    
    for (pe in seq(1:length(rownames(saved)))) {
      col_names <- colnames(annotate_result5[[dat]]@identification.result[[pe]])
      
      for (c in seq(1:nrow(annotate_result5[[dat]]@identification.result[[pe]]))) {
        line <- annotate_result5[[dat]]@identification.result[[pe]][c, ]
        x <- names(annotate_result5[[dat]]@identification.result[pe])
        x <- saved$MS2.spectra.name %in% x
        peak_number <- saved$MS1.peak.name[x]
        t <- unlist(c(peak_number, line))
        y <- rbind(y, t)
      }
    }
    
    if (iter == 1) {
      unione <- as.data.frame(y)
      colnames(unione)[1] <- "peak_number"
      unione_all <- rbind(unione_all, unione)
    } else {
      unione2 <- as.data.frame(y)
      colnames(unione2)[1] <- "peak_number"
      unione <- rbind(unione2[which(!unione2$peak_number %in% unione$peak_number), ], unione)
      unione_all <- rbind(unione_all, unione2)
    }
  }
  
  unione_all <- unione_all[!duplicated(unione_all$Compound.name) & !duplicated(unione_all$Lab.ID), ]
  unione_all$DATABASE <- "NA"
  unione_all$DATABASE[grep("*HMDB*", unione_all$Lab.ID)] <- "HMDB"
  unione_all$DATABASE[grep("MassBank*", unione_all$Lab.ID)] <- "MassBank"
  unione_all$DATABASE[grep("MONA*", unione_all$Lab.ID)] <- "MONA"
  
  unione <- unione[!duplicated(unione$Compound.name) & !duplicated(unione$Lab.ID), ]
  annot <- unione
  unione$DATABASE <- "FINAL_SELECTED"
  unione_all <- rbind(unione_all, unione)
  
  numbe <- as.numeric(gsub("peak", "", unione_all$peak_number))
  unione_all <- unione_all[order(numbe), ]
  
  write.table(unione_all,
              paste(directory2, "/", Cell_type, "_", "MS2_results.tsv", sep = ""),
              quote = FALSE, row.names = FALSE, sep = "\t")
  
  # Make variables globally available
  assign("directory2", directory2, envir = .GlobalEnv)
  assign("unione", unione, envir = .GlobalEnv)
  assign("unione_all", unione_all, envir = .GlobalEnv)
  assign("annot", annot, envir = .GlobalEnv)
}








# ---- Function to Generate MS2 Spectra Plots ----
generate_ms2_spectra_plots <- function(directory, directory2, Cell_type, args, MS2_databases, annotate_result5) {
  
  # Create subfolder for MS2 spectra
  dir.create(path = paste(directory2, "/MS2_SPECTRA", sep = ""), 
             showWarnings = TRUE, recursive = TRUE, mode = "0777")
  setwd(paste(directory2, "/MS2_SPECTRA", sep = ""))
  
  if (sum(MS2_databases == "HMDB") == 1) {
    Y <- grep("HMDB*", names(annotate_result5), ignore.case = TRUE)
    HMDA <- load(paste(directory, "/Integration/x_BiomiX_DATABASE/hmdb_database0.0.3.rda", sep = ""))
    peak_to_show <- names(annotate_result5[[Y]]@identification.result)
    matching <- annotate_result5[[Y]]@match.result
    show <- matching[matching$MS2.spectra.name %in% peak_to_show, ]
    ms2.plot1 <- ms2plot(object = annotate_result5[[Y]],
                         database = hmdb_database0.0.3,
                         which.peak = "all")
    file.rename("ms2_match_plot", "ms2_match_plot_HMDB")
  }
  
  if (sum(MS2_databases == "MASSBANK") == 1) {
    Y <- grep("MASSBANK*", names(annotate_result5), ignore.case = TRUE)
    HMDA <- load(paste(directory, "/Integration/x_BiomiX_DATABASE/massbank_database0.0.3.rda", sep = ""))
    peak_to_show <- names(annotate_result5[[Y]]@identification.result)
    matching <- annotate_result5[[Y]]@match.result
    show <- matching[matching$MS2.spectra.name %in% peak_to_show, ]
    ms2.plot1 <- ms2plot(object = annotate_result5[[Y]],
                         database = massbank_database0.0.3,
                         which.peak = "all")
    file.rename("ms2_match_plot", "ms2_match_plot_MASSBANK")
  }
  
  if (sum(MS2_databases == "MONA") == 1) {
    Y <- grep("MONA*", names(annotate_result5), ignore.case = TRUE)
    HMDA <- load(paste(directory, "/Integration/x_BiomiX_DATABASE/mona_database0.0.3.rda", sep = ""))
    peak_to_show <- names(annotate_result5[[Y]]@identification.result)
    matching <- annotate_result5[[Y]]@match.result
    show <- matching[matching$MS2.spectra.name %in% peak_to_show, ]
    ms2.plot1 <- ms2plot(object = annotate_result5[[Y]],
                         database = mona_database0.0.3,
                         which.peak = "all")
    file.rename("ms2_match_plot", "ms2_match_plot_MONA")
  }
  
  # Recreate output directory and reset working directory
  dir.create(path = paste(directory, "/Metabolomics/OUTPUT/", sep = ""), 
             showWarnings = TRUE, recursive = TRUE, mode = "0777")
  dir.create(path = paste(directory, "/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep = ""), 
             showWarnings = TRUE, recursive = TRUE, mode = "0777")
  directory2 <- paste(directory, "/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep = "")
  
  setwd(directory2)
  
  # Ensure important variables remain globally available
  assign("directory2", directory2, envir = .GlobalEnv)
  assign("peak_to_show", peak_to_show, envir = .GlobalEnv)
  assign("matching", matching, envir = .GlobalEnv)
  assign("show", show, envir = .GlobalEnv)
  assign("ms2.plot1", ms2.plot1, envir = .GlobalEnv)
}






# ---- Function to Retrieve or Generate Advanced Batch Annotation ----
retrieve_advanced_batch_annotation <- function(directory2, Cell_type, args, databases, masses_mode, Ion_mode, ducts, tolerance, MASS, RT) {
  
  # Check for existing annotation file
  if (file.exists(paste("Annotation_metabolites_", Cell_type, "_", args[1], "_vs_", args[2], ".tsv", sep = "")) == TRUE) {
    print("File available locally, using the local version")
    advanced_batch_df <- vroom(
      paste(directory2, "/Annotation_metabolites_", Cell_type, "_", args[1], "_vs_", args[2], ".tsv", sep = ""),
      delim = "\t", col_names = TRUE
    )
  } else {
    advanced_batch_df <- data.frame()
  }
  
  if (nrow(advanced_batch_df) == 0 || ncol(advanced_batch_df) == 0) {
    print("File is empty — using CMMR and CEUMASS mediator server to retrieve annotation.")
    
    
    #Allow to user to change the option based on their request  
    
    #  '["hmdb","metlin", "kegg", "lipidmaps"]'
    # '["M+H","M+Na","M+NH4","M+H-H2O"]'
    # databases <-list("hmdb", "metlin", "kegg", "lipidmaps" )
    # ducts <-list("M+H", "M+Na", "M+NH4", "M+H-H2O" )
    # masses_mode = "mz"
    # Ion_mode = "positive"
    # tolerance = 10
    # MASSI <- c(399.3367, 421.31686, 315.2424, 337.2234, 280.2402)
    # RT <- c(18.842525, 18.842525, 8.144917, 8.144917, 28.269503, 4.021555)
    
    
    # advanced_batch_df <- advanced_batch_search(
    #         cmm_url             = paste0(
    #                 'http://ceumass.eps.uspceu.es/mediator/api/v3/',
    #                 'advancedbatch'),
    #         chemical_alphabet   = 'all',
    #         modifiers_type      = 'none',
    #         metabolites_type    = 'all-except-peptides',
    #         databases           = databases,
    #         masses_mode         = masses_mode,
    #         ion_mode            = Ion_mode,
    #         adducts             = ducts,
    #         deuterium           = 'false',
    #         tolerance           = tolerance,
    #         tolerance_mode      = "ppm",
    #         #masses              = MASS,
    #         masses              = "[114.0907,245.0765]",
    #         all_masses          = '[]',
    #         #retention_times     = RT,
    #         retention_times     = "[2.316083333,0.948808333]",
    #         all_retention_times = '[]'
    # )
    
    advanced_batch_df <- advanced_batch_search(
      cmm_url = "https://ceumass.eps.uspceu.es/api/v3/advancedbatch",
      chemical_alphabet = "all",
      modifiers_type = "none",
      metabolites_type = "all-except-peptides",
      databases = databases,  #create a list
      masses_mode = as.character(masses_mode),
      ion_mode = as.character(Ion_mode),
      adducts = ducts,
      deuterium = FALSE,
      tolerance = as.numeric(tolerance),
      tolerance_mode = "ppm",
      masses = MASS,
      all_masses = list(),
      retention_times = RT ,
      all_retention_times = list(),
      composite_spectra = list(),
      all_composite_spectra = list()
    )
    
    head(advanced_batch_df)
    str(advanced_batch_df)
    
    print("File unavailable locally, generating the annotation from Ceu Mass Mediator database")
    write.table(
      advanced_batch_df,
      paste(directory2, "/Annotation_metabolites_", Cell_type, "_", args[1], "_vs_", args[2], ".tsv", sep = ""),
      quote = FALSE, row.names = FALSE, sep = "\t"
    )
    
  }
  
  # Make sure advanced_batch_df is available globally for downstream processing
  assign("advanced_batch_df", advanced_batch_df, envir = .GlobalEnv)
}




filter_by_biospecimen <- function(label, tat, metabolite_df, biospecimen_name) {
  # Check if this biospecimen should be used
  if (str_detect(label, fixed(biospecimen_name, ignore_case = TRUE))) {
    # Eliminate peaks not matching the biospecimen metabolite list
    to_eliminate2 <- NULL
    for (pe in unique(tat$name)) {
      sat <- which(tat$name == pe)
      if (sum(tat$HMDB[sat] %in% metabolite_df$HMDB_ID) == 0) {
        tat[sat, 9:ncol(tat)] <- NA
        to_eliminate2 <- append(to_eliminate2, sat[-1])
      }
    }
    tat <- tat[-to_eliminate2, ]
    
    # Merge with metabolite reference and filter out missing SMILES
    total <- merge(tat, metabolite_df, by.x = "HMDB", by.y = "HMDB_ID", all.x = TRUE)
    exit <- !is.na(total$identifier) & is.na(total$SMILES)
    total <- total[!exit, ]
    gc()
    return(list(total = total, filtering = "YES"))
  } else {
    return(NULL)  # No match, let caller handle fallback
  }
}






replace_ms1_with_ms2 <- function(annot, total, to_eliminate2 = NULL) {
  # This function replaces MS1 annotation with MS2 annotation when available.
  # It modifies total and to_eliminate2 exactly like the original code.
  
  for (yy in 1:length(unique(annot$peak_number))) {
    t <- annot$peak_number %in% unique(annot$peak_number)[yy]
    mix <- paste(annot$Compound.name[t], collapse = "/")
    mix2 <- paste(annot$Adduct[t], collapse = "/")
    mix3 <- paste(annot$mz.error[t], collapse = "/")
    
    if (sum(is.na(annot$HMDB.ID[t])) == 0) {
      x <- substr(total$HMDB, 7, 11) %in% substr(annot$HMDB.ID[t], 5, 9)
    } else {
      x <- FALSE
    }
    
    if (sum(x) == 0) {
      print(paste(unique(annot$peak_number)[yy], "not_annotated", sep = "_"))
      z <- which(total$NAME.x %in% unique(annot$peak_number)[yy])
      total[z, c(11:(ncol(total) - 2))] <- "NA"
      total[z[1], "name"] <- mix
      total[z[1], "adduct"] <- mix2
      total[z[1], "error_ppm"] <- mix3
      total[z[1], "MS2_annot"] <- mix
      total[z[1], "annot_level"] <- "Level_2"
      to_eliminate2 <- append(to_eliminate2, z[-1])
      
      # MS2 replaces MS1, keep only one annotation per peak.
    } else {
      t <- total$HMDB[x]
      t <- which(total$HMDB %in% t & total$NAME.x == unique(annot$peak_number)[yy])
      
      for (c in t) {
        print(paste(total$NAME.x[c], "already_annotated", sep = "_"))
        total$MS2_annot[c] <- mix
        total$name[c] <- mix
        total$adduct[c] <- mix2
        total$error_ppm[c] <- mix3
        total$annot_level[c] <- "Level_2"
      }
      
      out <- which(total$NAME.x == unique(annot$peak_number)[yy] &
                     total$annot_level == "Level_3")
      to_eliminate2 <- append(to_eliminate2, out)
      
      # MS2 annotated metabolites replace MS1 annotated ones.
    }
  }
  
  if (length(to_eliminate2) != 0) {
    total <- total[-to_eliminate2, ]
  }
  
  total$annot_level[is.na(total$name)] <- "Level_4"
  
  return(list(total = total, to_eliminate2 = to_eliminate2))
}




perform_dge_analysis_MS1_MS2 <- function(total, padju, LogFC, Cell_type, args, Metadata, directory2) {
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  total$padj <- as.numeric(total$padj)
  total$log2FC <- as.numeric(total$log2FC)
  
  # Get up- and down-regulated metabolites
  ALTI <- subset(total, padj < padju & log2FC > LogFC)
  ALTI <- ALTI[order(ALTI$padj), ]
  
  BASSI <- subset(total, padj < padju & log2FC < -LogFC)
  BASSI <- BASSI[order(BASSI$padj), ]
  
  x <- total$NAME.x %in% ALTI$NAME.x
  NO <- total[!x, ]
  x <- NO$NAME.x %in% BASSI$NAME.x
  NO <- NO[!x, ]
  
  # Write results to files
  write.table(ALTI$Name, file = paste(Cell_type, "_", args[1], "_", args[2], "_UP.tsv", sep = ""), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(BASSI$Name, file = paste(Cell_type, "_", args[1], "_", args[2], "_DOWN.tsv", sep = ""), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  write.table(ALTI, file = paste(Cell_type, "_", args[1], "_", args[2], "_METABUP.tsv", sep = ""), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(BASSI, file = paste(Cell_type, "_", args[1], "_", args[2], "_METABDOWN.tsv", sep = ""), 
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE), c("ID", "CONDITION")],
              file = paste(args[1], "_", args[2], "_", Cell_type, "_metadata.tsv", sep = ""),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  ALTI <- ALTI %>% distinct(NAME.x, .keep_all = TRUE)
  BASSI <- BASSI %>% distinct(NAME.x, .keep_all = TRUE)
  NO <- NO %>% distinct(NAME.x, .keep_all = TRUE)
  
  SKIP_ALTI <- nrow(ALTI) > 0
  SKIP_BASSI <- nrow(BASSI) > 0
  
  # Volcano plot
  p <- ggplot() +
    ggtitle(paste("DGE", Cell_type, "_", args[1], "_vs_", args[2], sep = "")) +
    theme(
      plot.title = element_text(color = "black", size = 14, face = "bold.italic", hjust = 0.5),
      axis.title.x = element_text(color = "black", size = 14, face = "bold"),
      axis.title.y = element_text(color = "black", size = 14, face = "bold")
    ) +
    ylab("-log10(padj)") + xlab("log2FC")
  
  if (SKIP_ALTI) {
    p <- p + geom_point(data = ALTI, aes(x = log2FC, y = -log10(padj)), color = "red") +
      geom_text_repel(data = ALTI[c(1:min(25, nrow(ALTI))), ],
                      aes(x = log2FC, y = -log10(padj), label = name),
                      hjust = 0, vjust = 0, size = 3, arrow = NULL, force = 10, force_pull = 1, max.overlaps = 10)
  }
  if (SKIP_BASSI) {
    p <- p + geom_point(data = BASSI, aes(x = log2FC, y = -log10(padj)), color = "blue") +
      geom_text_repel(data = BASSI[c(1:min(25, nrow(BASSI))), ],
                      aes(x = log2FC, y = -log10(padj), label = name),
                      hjust = 0, vjust = 0, size = 3, force = 20, force_pull = 1, max.overlaps = 10)
  }
  
  p <- p + geom_point(data = NO, aes(x = log2FC, y = -log10(padj)), color = "black")
  
  p2 <- p + geom_vline(xintercept = c(-LogFC, LogFC), col = "red") +
    geom_hline(yintercept = -log10(padju), col = "red")
  
  print(p2)
  dev.off()
  
  return(list(ALTI = ALTI, BASSI = BASSI))
}


generate_heatmap_MS1_MS2 <- function(ALTI, BASSI, matrix, Heatmap_genes, Metadata, Cell_type, args, directory2) {
  if (length(c(ALTI$NAME.x, BASSI$NAME.x)) > 2) {
    if (nrow(ALTI) > Heatmap_genes) {
      ALTI <- ALTI[1:Heatmap_genes, ]
    }
    if (nrow(BASSI) > Heatmap_genes) {
      BASSI <- BASSI[1:Heatmap_genes, ]
    }
    
    Heat <- matrix[, colnames(matrix) %in% c(ALTI$NAME.x, BASSI$NAME.x), drop = FALSE]
    samples <- rownames(Heat)
    Heat <- as.data.frame(Heat)
    Heat <- apply(Heat, 2, as.numeric)
    rownames(Heat) <- samples
    Heat <- t(Heat)
    Heat[Heat == 0] <- 1
    Heat <- apply(Heat, 2, log10)
    
    annot_vis <- rownames(Heat)
    pool <- rbind(ALTI, BASSI)
    
    SAVED <- NULL
    for (ix in annot_vis) {
      iu <- pool$NAME.x %in% ix
      Gene <- as.character(pool$Name[which(iu)])
      SAVED <- append(Gene, SAVED)
    }
    SAVED <- rev(as.character(SAVED))
    
    for (ix in 1:length(samples)) {
      if (!is.na(SAVED[ix])) {
        rownames(Heat)[ix] <- SAVED[ix]
      }
    }
    
    library(ComplexHeatmap)
    library(circlize)
    
    setwd(directory2)
    pdf(file = paste("Heatmap_top_genes_", Cell_type, "_", args[2], "_vs_", args[1], ".pdf", sep = ""))
    
    col_fun <- colorRamp2(c(min(Heat, na.rm = TRUE),
                            mean(colMeans(Heat, na.rm = TRUE)),
                            max(Heat, na.rm = TRUE)),
                          c("blue", "black", "yellow"))
    
    t <- c("CTRL" = "blue", "SLE" = "red")
    attr(t, "names")[1] <- args[2]
    attr(t, "names")[2] <- args[1]
    
    ha <- HeatmapAnnotation(condition = Metadata$CONDITION,
                            col = list(condition = t))
    
    Heat[is.na(Heat)] <- 0
    
    p <- Heatmap(Heat, km = 2, name = "SD_score", col = col_fun,
                 clustering_distance_rows = "pearson", clustering_method_rows = "complete",
                 clustering_method_columns = "ward.D",
                 row_dend_width = unit(0.5, "cm"), column_dend_height = unit(60, "mm"),
                 column_names_gp = grid::gpar(fontsize = 6),
                 row_names_gp = grid::gpar(fontsize = 8),
                 top_annotation = ha)
    
    print(p)
    
    colnames(Heat) <- Metadata$CONDITION
    ha <- HeatmapAnnotation(condition = Metadata$CONDITION,
                            col = list(condition = t))
    
    p <- Heatmap(Heat, km = 2, name = "SD_score", col = col_fun,
                 clustering_distance_rows = "pearson", clustering_method_rows = "complete",
                 clustering_method_columns = "ward.D",
                 row_dend_width = unit(0.5, "cm"), column_dend_height = unit(60, "mm"),
                 column_names_gp = grid::gpar(fontsize = 6),
                 row_names_gp = grid::gpar(fontsize = 8),
                 top_annotation = ha)
    
    print(p)
    dev.off()
  }
}




run_metpath_pipeline_MS1_MS2 <- function(total,
                                     saved,
                                     directory2,
                                     Cell_type,
                                     args,
                                     hmdb_pathway,
                                     kegg_hsa_pathway,
                                     COMMAND_ADVANCED) {
  
  # Filter significant hits (still using p_val as in your code)
  query_id <- total[which(total$p_val < 0.05), ]  # TODO: consider replacing with p.adj if needed
  query_id_saved <- saved[which(saved$p_val < 0.05), ]
  assign("query_id", query_id, envir = .GlobalEnv)
  assign("query_id_saved", query_id_saved, envir = .GlobalEnv)
  
  # Run pathway analysis only if significant metabolites exist
  if (nrow(query_id) != 0) {
    
    ## --- HMDB Pathway Analysis ---
    dir.create(path = paste(directory2, "/Pathway_analysis/", "HMDB_", Cell_type, "_", args[1], "_vs_", args[2], sep = ""),
               showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2, "/Pathway_analysis/", "HMDB_", Cell_type, "_", args[1], "_vs_", args[2], sep = ""))
    
    pathway_class_HMDB <- metpath::pathway_class(hmdb_pathway)
    pathway_class_KEGG <- metpath::pathway_class(kegg_hsa_pathway)
    
    pdf(file = paste("Pathway_analysis_HMDB_", args[1], "_", args[2], "_", Cell_type, ".pdf", sep = ""))
    
    gc()
    remain_idx <- which(unlist(pathway_class_HMDB) == "Metabolic;primary_pathway")
    hmdb_pathway <- hmdb_pathway[remain_idx]
    
    result <- enrich_hmdb(query_id = unique(query_id$HMDB),
                          query_type = "compound",
                          id_type = "HMDB",
                          pathway_database = hmdb_pathway,
                          only_primary_pathway = TRUE,
                          p_cutoff = 0.05,
                          p_adjust_method = "BH",
                          threads = as.numeric(COMMAND_ADVANCED[3, 3]))
    
    if (length(result) != 0) {
      x <- enrich_bar_plot(object = result,
                           x_axis = "p_value_adjust",
                           cutoff = 1.1,
                           top = 10)
      print(x)
      
      x <- enrich_scatter_plot(object = result)
      print(x)
      
      write.table(result@result,
                  paste(directory2, "/Pathway_analysis/", "HMDB_", Cell_type, "_", args[1], "_vs_", args[2], "/HMDB_table_results.tsv", sep = ""),
                  quote = FALSE, row.names = FALSE, sep = "\t")
      dev.off()
    }
    
    ## --- KEGG Pathway Analysis ---
    gc()
    dir.create(path = paste(directory2, "/Pathway_analysis/", "KEGG_", Cell_type, "_", args[1], "_vs_", args[2], sep = ""),
               showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2, "/Pathway_analysis/", "KEGG_", Cell_type, "_", args[1], "_vs_", args[2], sep = ""))
    
    pdf(file = paste("Pathway_analysis_KEGG_", args[1], "_", args[2], "_", Cell_type, ".pdf", sep = ""))
    
    head(pathway_class_KEGG)
    remain_idx <- pathway_class_KEGG %>%
      unlist() %>%
      stringr::str_detect("Disease") %>%
      `!`() %>%
      which()
    
    pathway_database <- kegg_hsa_pathway[remain_idx]
    
    result <- enrich_kegg(query_id = unique(query_id$KEGG),
                          query_type = "compound",
                          id_type = "KEGG",
                          pathway_database = pathway_database,
                          p_cutoff = 0.05,
                          p_adjust_method = "BH",
                          threads = as.numeric(COMMAND_ADVANCED[3, 3]))
    
    if (length(result) != 0) {
      x <- enrich_bar_plot(object = result,
                           x_axis = "p_value_adjust",
                           cutoff = 1.1,
                           top = 10)
      print(x)
      
      x <- enrich_scatter_plot(object = result)
      print(x)
      
      write.table(result@result,
                  paste(directory2, "/Pathway_analysis/", "KEGG_", Cell_type, "_", args[1], "_vs_", args[2], "/KEGG_table_results.tsv", sep = ""),
                  quote = FALSE, row.names = FALSE, sep = "\t")
      dev.off()
    }
    
  } else {
    print("NO STATISTICAL SIGNIFICANT PEAK FOR METABOLITE PATHWAY ANALYSIS")
  }
}

  
  
  
  run_metaboanalyst_pipeline_MS1_MS2 <- function(
    query_id,
    query_id_saved,
    total,
    Cell_type,
    args,
    directory,
    directory2
  ) {
    # Ensure results are written in the global environment if needed
    # (no explicit return required since function creates files, not new variables)
    
    #Enrichment analysis
    dir.create(path = paste(directory2,"/MetaboAnalyst/", "Enrichment_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/", "Enrichment_Analysis", sep =""))
    
    numbe<-as.numeric(gsub("peak", "", query_id$NAME.x))
    query_id<-query_id[order(numbe),]
    
    write.table(x= query_id$HMDB[!is.na(query_id$HMDB)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv", sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id$KEGG[!is.na(query_id$KEGG )] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id$name[!is.na(query_id$name )] , file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    
    numbe<-as.numeric(gsub("peak", "", query_id_saved$NAME.x))
    query_id_saved<-query_id_saved[order(numbe),]
    write.table(x= query_id_saved[,c(3,1,21,14)] , file= paste("HMDB_KEGG_Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],"_complete",".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
    
    
    #Pathway analysis
    dir.create(path = paste(directory2,"/MetaboAnalyst/", "Pathway_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/", "Pathway_Analysis", sep =""))
    
    write.table(x= query_id$HMDB[!is.na(query_id$HMDB)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id$KEGG[!is.na(query_id$KEGG )] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id$name[!is.na(query_id$name )] , file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= ",", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id_saved[,c(3,1,21,14)] , file= paste("HMDB_KEGG_Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],"_complete.tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
    
    
    #Joint-Pathway analysis
    dir.create(path = paste(directory2,"/MetaboAnalyst/", "Joint_Pathway_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/", "Joint_Pathway_Analysis", sep =""))
    
    query_id_select <- query_id[,c("HMDB", "KEGG", "name", "log2FC")]
    
    write.table(x= query_id_select[,c(1,4)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id_select[,c(2,4)] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id_select[,c(3,4)], file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    
    #LATE INTEGRATION TRANSCRIPTOMICS
    regeX <- paste("*",args[1],"_vs_",args[2], sep="")
    regeX2 <- paste(args[1],"_vs_",args[2], sep="")
    files <- grep(regeX,list.files(paste(directory,"/Transcriptomics/OUTPUT",sep="")),value=TRUE)
    
    if (length(files) != 0){
      for (fil in 1:length(files)){
        print(fil)
        files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2, sep="")),value=TRUE)
        if (length(files2) == 0) {
          message("No matching trascriptomics files found or not enough files — skipping to next")
          next
        }
        LIST_GENE<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[1], sep=""), delim="\t")
        LIST_GENE2<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[2], sep=""), delim="\t")
        LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
        
        Sources<-str_split(files2[1], "_")[[1]][1]
        print(Sources)
        
        dir.create(path = paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Transcriptomes_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Transcriptomes_availables/", Sources, sep =""))
        
        write.table(x= LIST_GENE3[c("Gene.name", "log2FoldChange")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
      }
    }
    
    #LATE INTEGRATION METHYLOMICS
    regeX <- paste("*",args[1],"_vs_",args[2], sep="")
    regeX2 <- paste(args[1],"_vs_",args[2], sep="")
    files <- grep(regeX,list.files(paste(directory,"/Methylomics/OUTPUT",sep="")),value=TRUE)
    
    if (length(files) != 0){
      for (fil in 1:length(files)){
        print(fil)
        files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Methylomics/OUTPUT/", files[fil], sep="")),value=TRUE)
        if (length(files2) == 0) {
          message("No matching methylation files found or not enough files — skipping to next")
          next
        }
        LIST_GENE<-vroom(paste(directory,"/Methylomics/OUTPUT/", files[fil],"/",files2[1], sep=""), delim="\t")
        LIST_GENE2<-vroom(paste(directory,"/Methylomics/OUTPUT/", files[fil],"/",files2[2], sep=""), delim="\t")
        LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
        
        Sources<-str_split(files2[1], "_")[[1]][1]
        print(Sources)
        
        dir.create(path = paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Methylomics_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Methylomics_availables/", Sources, sep =""))
        
        write.table(x= LIST_GENE3[c("gene", "logFC")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
      }
    }
    
    #Network analysis
    dir.create(path = paste(directory2,"/MetaboAnalyst/", "Network_Analysis", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/", "Network_Analysis", sep =""))
    
    query_id_select <- query_id[,c("HMDB", "KEGG", "name", "log2FC")]
    
    write.table(x= query_id_select[,c(1,4)] , file= paste("HMDB_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id_select[,c(2,4)] , file= paste("KEGG_ID_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    write.table(x= query_id_select[,c(3,4)], file= paste("Compound_names_",Cell_type,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
    
    #LATE INTEGRATION TRANSCRIPTOMICS
    regeX <- paste("*",args[1],"_vs_",args[2], sep="")
    regeX2 <- paste(args[1],"_vs_",args[2], sep="")
    files <- grep(regeX,list.files(paste(directory,"/Transcriptomics/OUTPUT",sep="")),value=TRUE)
    
    if (length(files) != 0){
      for (fil in 1:length(files)){
        print(fil)
        files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2, sep="")),value=TRUE)
        if (length(files2) == 0) {
          message("No matching transcriptomics files found or not enough files — skipping to next")
          next
        }
        LIST_GENE<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[1], sep=""), delim="\t")
        LIST_GENE2<-vroom(paste(directory,"/Transcriptomics/OUTPUT/", files[fil],"/",regeX2,"/",files2[2], sep=""), delim="\t")
        LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
        
        Sources<-str_split(files2[1], "_")[[1]][1]
        print(Sources)
        
        dir.create(path = paste(directory2,"/MetaboAnalyst/Network_Analysis/Transcriptomes_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/Network_Analysis/Transcriptomes_availables/", Sources, sep =""))
        
        write.table(x= LIST_GENE3[c("Gene.name", "log2FoldChange")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
      }
    }
    
    #LATE INTEGRATION METHYLOMICS
    regeX <- paste("*",args[1],"_vs_",args[2], sep="")
    regeX2 <- paste(args[1],"_vs_",args[2], sep="")
    files <- grep(regeX,list.files(paste(directory,"/Methylomics/OUTPUT",sep="")),value=TRUE)
    
    if (length(files) != 0){
      for (fil in 1:length(files)){
        print(fil)
        files2 <- grep("GENES[A-Z]*",list.files(paste(directory,"/Methylomics/OUTPUT/", files[fil], sep="")),value=TRUE)
        if (length(files2) == 0) {
          message("No matching methylomics files found or not enough files — skipping to next")
          next
        }
        LIST_GENE<-vroom(paste(directory,"/Methylomics/OUTPUT/", files[fil],"/",files2[1], sep=""), delim="\t")
        LIST_GENE2<-vroom(paste(directory,"/Methylomics/OUTPUT/", files[fil],"/",files2[2], sep=""), delim="\t")
        LIST_GENE3<-rbind(LIST_GENE, LIST_GENE2)
        
        Sources<-str_split(files2[1], "_")[[1]][1]
        print(Sources)
        
        dir.create(path = paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Methylomics_availables/", Sources, sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/MetaboAnalyst/Joint_Pathway_Analysis/Methylomics_availables/", Sources, sep =""))
        
        write.table(x= LIST_GENE3[c("gene", "logFC")], file= paste("GENE_ID_",Sources,"_",args[1],"_vs_",args[2],".tsv",sep="")  ,sep= "\t", row.names = FALSE, col.names = FALSE,  quote = FALSE)
      }
    }
    
    #MSEA ANALYSIS
    dir.create(path = paste(directory2,"/MetaboAnalyst/", "Functional_analysis_MSEA", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
    setwd(paste(directory2,"/MetaboAnalyst/", "Functional_analysis_MSEA", sep =""))
    
    MSEA <- total[,c("m/z","p_val", "log2FC", "RT_min" )]
    colnames(MSEA) <- c("m.z","p.value", "t.score", "r.t" )
    
    write.table(x= MSEA, file= paste("MSEA_",args[1],"_vs_",args[2],".txt",sep="")  ,sep= "\t", row.names = FALSE, col.names = TRUE,  quote = FALSE)
    gc()
  }
  


