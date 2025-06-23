
# Check type of gene name used (ENSEMBL or GENE SYMBOL)

library(tidyverse)

process_gene_annotation <- function(Matrix, genes, directory) {
        
        # Load mart file
        mart_path <- file.path(directory, "Integration", "x_BiomiX_DATABASE", "mart_export_37.txt")
        Mart <- read.table(mart_path, sep = ",", header = TRUE)
        
        Matrix <- as.data.frame(Matrix)
        Matrix$X <- genes
        
        if (str_detect(rownames(Matrix)[1], "ENSG")) {
                # ENSG-based annotation
                DGE <- merge(Matrix, Mart, by.x = "X", by.y = "Gene.stable.ID")
                genes_name <- DGE$X
                genes_name_C <- DGE$Gene.name
                DGE2 <- DGE %>%
                        column_to_rownames(var = "X")
                DGE2 <- DGE2[, -ncol(DGE2)]
                GENE_ANNOTATION <- "ENSEMBL"
                
        } else {
                # Gene name-based annotation
                DGE <- Matrix
                DGE$Gene.name <- DGE$X
                genes_name <- DGE$X
                genes_name_C <- DGE$Gene.name
                
                DGE2 <- DGE[, -(ncol(DGE)-1:0)] %>% 
                        as.matrix()
                
                rownames(DGE2) <- genes_name
                DGE2 <- as.data.frame(DGE2)
                GENE_ANNOTATION <- "GENE_NAME"
        }
        
        return(list(
                DGE2 = DGE2,
                genes_name = genes_name,
                genes_name_C = genes_name_C,
                GENE_ANNOTATION = GENE_ANNOTATION
        ))
}


#FILTERING SAMPLES BASED ON METADATA CRITERIA SELECTED. 

library(tidyverse)

filter_metadata_and_expression <- function(Metadata_Bcell, DGE2, COMMAND_ADVANCED, directory) {
        
        # Set working directory (optional, for backward compatibility)
        setwd(file.path(directory, "Transcriptomics", "INPUT"))
        
        # Identify FILTERING columns
        meta_filter_indices <- grep(".*FILTERING.*", colnames(COMMAND_ADVANCED))
        
        for (i in seq_along(meta_filter_indices)) {
                idx <- meta_filter_indices[i]
                filter_value <- COMMAND_ADVANCED[3, idx]
                
                if (!is.na(filter_value)) {
                        colname <- as.character(COMMAND_ADVANCED[1, idx])
                        filter_type <- as.character(COMMAND_ADVANCED[2, idx])
                        condition <- as.character(filter_value)
                        
                        if (filter_type == "Numerical") {
                                to_filter <- as.numeric(unlist(Metadata_Bcell[[colname]]))
                                
                                # Extract operator and numeric threshold
                                symbol <- str_extract(condition, "[><=]+")
                                threshold <- as.numeric(str_remove_all(condition, "[><= ]"))
                                
                                # Define comparison function
                                comparison_operator <- switch(symbol,
                                                              "<"  = function(a, b) a < b,
                                                              ">"  = function(a, b) a > b,
                                                              "==" = function(a, b) a == b,
                                                              ">=" = function(a, b) a >= b,
                                                              "<=" = function(a, b) a <= b,
                                                              stop("Invalid operator"))
                                
                                keep <- comparison_operator(to_filter, threshold)
                                Metadata_Bcell <- Metadata_Bcell[keep, ]
                                DGE2 <- DGE2[, keep, drop = FALSE]
                                
                        } else if (filter_type == "Factors") {
                                to_filter <- as.character(unlist(Metadata_Bcell[[colname]]))
                                
                                # Extract operator and factor value
                                symbol <- str_extract(condition, "==|!=")
                                threshold <- str_remove_all(condition, "==|!=| ")
                                
                                comparison_operator <- switch(symbol,
                                                              "==" = function(a, b) a == b,
                                                              "!=" = function(a, b) a != b,
                                                              stop("Invalid operator"))
                                
                                keep <- comparison_operator(to_filter, threshold)
                                Metadata_Bcell <- Metadata_Bcell[keep, ]
                                DGE2 <- DGE2[, keep, drop = FALSE]
                        }
                }
        }
        
        return(list(
                Metadata_Bcell = Metadata_Bcell,
                DGE2 = DGE2
        ))
}

#Preview function

launch_qc_preview <- function(DGE2, Metadata_Bcell, directory, browser_analysis) {
        # Set up numeric matrix for visualization
        Samples_preview <- colnames(DGE2)
        numeric_data <- t(DGE2)
        numeric_data <- apply(numeric_data, 2, as.numeric)
        rownames(numeric_data) <- Samples_preview
        numeric_data <- as.data.frame(numeric_data)
        
        # Prepare metadata
        metadata <- Metadata_Bcell
        colnames(metadata)[colnames(metadata) == "CONDITION"] <- "condition"
        
        # Load Shiny app
        source(file.path(directory, "BiomiX_preview.r"))
        browser_path <- readLines(file.path(directory, "_INSTALL", "CHOISE_BROWSER_pre-view"), n = 1)
        
        options(browser = browser_path)
        message("Pre QC data visualization. If the browser does not open automatically, copy and paste the link into a browser.")
        message("Once the QC window opens, apply your filters and click 'close app' to continue.")
        message("NOTE: If the browser does not open, set up the correct path in CHOISE_BROWSER_pre-view.")
        
        # Run the Shiny App
        Preview <- shiny::runApp(runShinyApp(numeric_data, metadata), launch.browser = TRUE)
        
        # Update global variables
        DGE2_updated <- as.data.frame(t(Preview$matrix))
        Metadata_Bcell_updated <- Preview$metadata
        colnames(Metadata_Bcell_updated)[colnames(Metadata_Bcell_updated) == "condition"] <- "CONDITION"
        
        return(list(DGE2 = DGE2_updated, Metadata_Bcell = Metadata_Bcell_updated))
}



#STATEMENT IF DATA ARE NORMALIZED OR NOT, IT LOOKS IF NUMERS ARE FLOATS OR INTEGER
#Normalization, process

library(DESeq2)
library(tidyverse)

is_float <- function(x) {
        is.numeric(x) && any(x %% 1 != 0)
}

handle_normalization_and_export <- function(DGE2, Metadata_Bcell, Cell_type, args, directory) {
        
        # Detect if normalized
        if (any(sapply(DGE2[1:10, 1:min(10, ncol(DGE2))], is_float))) {
                NORMALIZATION <- "YES"
                message(paste("Automatic detection of normalized data. Are the data normalized? ", NORMALIZATION))
                
                # Prepare matrix
                DGE3 <- DGE2
                
                DGE2 <- as.data.frame(DGE2) %>%
                        rownames_to_column(var = "ID")

                # Create output directory
                out_dir <- file.path(directory, "Integration", "INPUT", paste0(Cell_type, "_", args[1], "_vs_", args[2]))
                dir.create(out_dir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
                
                # Export metadata
                write_tsv(Metadata_Bcell, file.path(out_dir, paste0("Metadata_", Cell_type, "_", args[1], ".tsv")))
                
                DGE3_out <- as.data.frame(DGE3) %>%
                        rownames_to_column(var = "Gene.name")
                
                # Rename columns to match metadata
                colnames(DGE3_out)<-  c("ID",Metadata_Bcell$ID)
                write_tsv(DGE3_out, file.path(out_dir, paste0(Cell_type, "_", args[1], "_vs_", args[2], "_normalized_vst.tsv")))
                
                return(list(NORMALIZATION = NORMALIZATION, DGE2 = DGE2, DGE3 = DGE3_out ))
        } else {
                NORMALIZATION <- "NO"
                message(paste("Automatic detection of normalized data. Are the data normalized? ", NORMALIZATION))
                
                # -------- Perform DESeq2 Normalization & VST --------
                DGE2 <- apply(DGE2, 2, round)
                
                dds <- DESeqDataSetFromMatrix(
                        countData = DGE2,
                        colData = Metadata_Bcell,
                        design = ~1
                )
                
                dds <- estimateSizeFactors(dds)
                
                normalized_counts <- counts(dds, normalized = TRUE)
                normalized_counts_vst <- vst(assay(dds))
                
                # Format output matrices
                DGE2_out <- as.data.frame(normalized_counts) %>%
                        rownames_to_column(var = "ID")
                
                DGE3_out <- as.data.frame(assay(dds)) %>%
                        rownames_to_column(var = "ID")
                
                DGE4_out <- as.data.frame(normalized_counts_vst) %>%
                        rownames_to_column(var = "ID")
                
                # Output directory
                out_dir <- file.path(directory, "Integration", "INPUT", paste0(Cell_type, "_", args[1], "_vs_", args[2]))
                dir.create(out_dir, recursive = TRUE, showWarnings = FALSE, mode = "0777")
                
                # Export metadata
                write_tsv(Metadata_Bcell, file.path(out_dir, paste0("Metadata_", Cell_type, "_", args[1], ".tsv")))
                
                # Rename expression matrix columns
                colnames(DGE2_out) <- c("ID",Metadata_Bcell$ID)
                colnames(DGE3_out) <- c("ID",Metadata_Bcell$ID)
                colnames(DGE4_out) <- c("ID",Metadata_Bcell$ID)
                
                # Export files
                write_tsv(DGE2_out, file.path(out_dir, paste0(Cell_type, "_", args[1], "_vs_", args[2], "_normalized.tsv")))
                write_tsv(DGE3_out, file.path(out_dir, paste0(Cell_type, "_", args[1], "_vs_", args[2], "_unnormalized.tsv")))
                write_tsv(DGE4_out, file.path(out_dir, paste0(Cell_type, "_", args[1], "_vs_", args[2], "_normalized_vst.tsv")))
                
                return(list(NORMALIZATION = NORMALIZATION, DGE2 = DGE2_out, DGE3 = DGE3_out, dds = dds))
        }
}




#### HEATMAP SECTION ----

generate_heatmap_signature <- function(DGE2, dds, Metadata_Bcell, Gene_panel, GENE_ANNOTATION, Mart, directory, Cell_type, args, COMMAND_ADVANCED) {
        library(ComplexHeatmap)
        library(circlize)
        library(readxl)
        library(vroom)
        library(dplyr)
        library(stringr)
        
        directory2 <- file.path(directory, "Transcriptomics/INPUT")
        setwd(directory2)
        
        # Load genes
        if (COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3] != "X") {
                genes <- if (grepl("\\.xlsx$|\\.xls$", Gene_panel)) {
                        read_excel(Gene_panel)$GENES_FOR_SUBPOPULATION
                } else {
                        vroom(Gene_panel, delim = "\t", col_names = TRUE)$GENES_FOR_SUBPOPULATION
                }
        } else {
                print("No panel selected")
                return(Metadata_Bcell)
        }
        
        directory2 <- file.path(directory, "Integration/INPUT", paste0(Cell_type, "_", args[1], "_vs_", args[2]))
        setwd(directory2)
        
        # Use appropriate normalization
        DGE <- if (NORMALIZATION == "NO") {
                counts(dds, normalized = TRUE)
        } else {
                DGE2 %>% column_to_rownames("ID")
        }
        
        DGE <- as.data.frame(DGE)
        DGE$X <- rownames(DGE)
        
        if (GENE_ANNOTATION == "GENE_NAME") {
                DGE$Gene.name <- DGE$X
                DGE <- DGE[, !grepl("X", colnames(DGE))]
        } else {
                DGE <- merge(DGE, Mart, by.x = "X", by.y = "Gene.stable.ID")
        }
        
        Panel <- which(DGE$Gene.name %in% genes)
        DGE <- DGE[Panel, ]
        rownames(DGE) <- DGE$Gene.name
        DGE <- DGE[, !colnames(DGE) %in% c("Gene.name", "X")]
        
        # Create z-score matrix
        IFN_score <- as.data.frame(matrix(0, ncol = ncol(DGE), nrow = nrow(DGE)))
        colnames(IFN_score) <- colnames(DGE)
        rownames(IFN_score) <- genes
        
        HC <- DGE[, Metadata_Bcell$CONDITION == args[2]]
        for (gene in 1:nrow(DGE)) {
                mean_val <- mean(as.numeric(HC[gene, ]))
                sd_val <- sd(as.numeric(HC[gene, ]))
                IFN_score[gene, ] <- (as.numeric(DGE[gene, ]) - mean_val) / sd_val
        }
        
        # Subset IFN_score for the relevant samples
        condition_filter <- Metadata_Bcell$CONDITION %in% c(args[1], args[2])
        IFN_score <- IFN_score[, condition_filter]
        
        # Automatic positivity detection
        positivity1 <- colSums(IFN_score > 2, na.rm = TRUE) >= as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING[1])
        positivity2 <- colSums(IFN_score > 1, na.rm = TRUE) >= as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING[2])
        Metadata_Bcell$condition <- ifelse(positivity1 | positivity2, "pos", "neg")
        
        # Heatmap coloring setup
        col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "black", "yellow"))
        ha <- HeatmapAnnotation(condition = Metadata_Bcell$condition[condition_filter],
                                col = list(condition = c("neg" = "blue", "pos" = "red")))
        
        pdf(file = paste0("Gene_panel_subgroups_", Cell_type, "_", args[1], "_vs_", args[2], ".pdf"))
        heatmap <- Heatmap(IFN_score,
                           km = 2,
                           name = "SD_score",
                           col = col_fun,
                           clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],
                           clustering_method_rows = "complete",
                           clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2],
                           row_dend_width = unit(0.5, "cm"),
                           column_dend_height = unit(60, "mm"),
                           column_names_gp = grid::gpar(fontsize = 6),
                           row_names_gp = grid::gpar(fontsize = 8),
                           top_annotation = ha)
        print(heatmap)
        dev.off()
        
        # Optional validation
        if ("MARKER" %in% colnames(Metadata_Bcell)) {
                directory2 <- file.path(directory, "Transcriptomics/INPUT")
                setwd(directory2)
                
                sink("Validation_MARKER_vs_GENE_PANEL.tsv")
                
                for (group in c(args[1], args[2])) {
                        group_data <- Metadata_Bcell %>% filter(CONDITION == group)
                        measured <- group_data %>% filter(MARKER > 0 | MARKER < 0)
                        cat(sprintf("Total %s samples: %d\n", group, nrow(group_data)))
                        cat(sprintf("%s with IFN measurement: %d\n", group, nrow(measured)))
                        
                        pos_classified <- measured %>% filter(condition == "pos")
                        neg_classified <- measured %>% filter(condition == "neg")
                        cat(sprintf("Positives correctly classified: %d / %d\n",
                                    nrow(pos_classified %>% filter(MARKER > 1.25)),
                                    nrow(group_data %>% filter(MARKER > 1.25))))
                        cat(sprintf("Negatives correctly classified: %d / %d\n",
                                    nrow(neg_classified %>% filter(MARKER < 1.25)),
                                    nrow(group_data %>% filter(MARKER < 1.25))))
                }
                
                sink()
        } else {
                print("No validation marker added")
        }
        
        directory2 <- file.path(directory, "Integration/INPUT", paste0(Cell_type, "_", args[1], "_vs_", args[2]))
        setwd(directory2)
        return(Metadata_Bcell)
}


#### SAVING NORMALIZED DATA ----

#Function to replace the negative with healthy signature
prepare_metadata_conditions <- function(Metadata) {
        num <- which(Metadata$CONDITION == args[2] & Metadata$condition == "neg")
        Metadata$condition[num] <- "healthy"
        return(Metadata)
}

#Function to remove the CTRL positive for the marker
filter_samples_by_condition <- function(DGE2, DGE3, Metadata) {
        num <- which(Metadata$CONDITION == args[2] & Metadata$condition == "pos")
        t<-as.character(names(num))
        if (length(num) != 0){
        DGE2 <- DGE2 %>% dplyr::select(!names(num))
        DGE3 <- DGE3 %>% dplyr::select(!names(num))
        Metadata <- Metadata[-num, ]
        }else{
                DGE2 <- DGE2
                DGE3 <- DGE3
                Metadata <- Metadata
        }
        return(list(DGE2 = DGE2, DGE3 = DGE3, Metadata = Metadata))
}

format_expression_matrix <- function(DGE2, Metadata, genes_name, normalized) {
        DGE2<- DGE2 %>% t() %>% as.data.frame() %>% dplyr::slice(-1)
        # if (normalized == "NO") {
        #         DGE2 <- DGE2[-1, ]
        # }
        DGE2$condition <- Metadata$condition
        colnames(DGE2) <- c(genes_name, "condition")
        rownames(DGE2) <- Metadata$ID
        DGE2$samples <- Metadata$ID
        return(DGE2)
}

save_outputs <- function(Metadata, DGE3, directory, cell_type, group1, group2) {
        setwd(directory)
        write.table(Metadata, file = paste("Metadata_", cell_type, "_", group1, "_vs_", group2, "_final.tsv", sep = ""),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        write.table(DGE3, file = paste("gene_count_", cell_type, "_", group1, "_vs_", group2, "_unnormalized.tsv", sep = ""),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}






#### GENERATION MATRIX FOR DATA INTEGRATION ----


generate_mofa_input <- function(directory, Cell_type, args, NORMALIZATION, GENE_ANNOTATION, COMMAND_ADVANCED) {
        library(tidyverse)
        library(vroom)
        library(ggplot2)
        library(dplyr)
        
        # Set working directory
        directory2 <- file.path(directory, "Integration", "INPUT", paste0(Cell_type, "_", args[1], "_vs_", args[2]))
        setwd(directory2)
        
        # Load normalized data
        if (NORMALIZATION == "NO") {
                print("PATH NO NORMALIZED")
                
                normalizzati <- vroom(paste0(Cell_type, "_", args[1], "_vs_", args[2], "_normalized_vst.tsv"), delim = "\t")
                normalizzati_size <- vroom(paste0(Cell_type, "_", args[1], "_vs_", args[2], "_normalized.tsv"), delim = "\t")
                
                # Set gene names as rownames
                genes <- normalizzati$ID
                normalizzati <- normalizzati %>% column_to_rownames("ID")
                normalizzati_size <- normalizzati_size %>% column_to_rownames("ID")
                
                # Filter out genes with low expression (expressed in <= half the samples)
                expressed_mask <- rowSums(normalizzati_size > 0) > (ncol(normalizzati_size) / 2)
                normalizzati <- normalizzati[expressed_mask, ]
                normalizzati_size <- normalizzati_size[expressed_mask, ]
                gene <- rownames(normalizzati)
                
        } else {
                print("PATH YES NORMALIZED")
                
                normalizzati <- vroom(paste0(Cell_type, "_", args[1], "_vs_", args[2], "_normalized_vst.tsv"), delim = "\t")
                genes <- normalizzati$ID
                normalizzati <- normalizzati  %>% column_to_rownames("ID")
                
                # Filter low expression genes
                expressed_mask <- rowSums(normalizzati > min(normalizzati)) > (ncol(normalizzati) / 2)
                normalizzati <- normalizzati[expressed_mask, ]
                gene <- rownames(normalizzati)
        }
        
        print("Genes before filtering:")
        print(length(genes))
        print("Genes after filtering:")
        print(nrow(normalizzati))
        
        #================ VARIANCE FILTERING ================
        
        normalizzati <- as.data.frame(lapply(normalizzati, as.numeric))
        varianze <- rowVars(as.matrix(normalizzati), na.rm = TRUE)
        # varianze <- apply(as.matrix(normalizzati), 1, var) #For some reason not working on Windows terminal

        rownames(normalizzati) <- gene
        # Add ID and variance back into the data frame
        normalizzati <- normalizzati %>%
                as.data.frame() %>%
                rownames_to_column("ID") %>%
                mutate(variance = log10(varianze))
        
        if (NORMALIZATION == "NO") {
                normalizzati_size <- normalizzati_size %>%
                        as.data.frame() %>%
                        rownames_to_column("ID") %>%
                        mutate(variance = log10(varianze))
        }
        
        # Merge with gene annotation if needed
        if (GENE_ANNOTATION == "GENE_NAME") {
                MART <- read.table(file.path(directory, "Integration", "x_BiomiX_DATABASE", "mart_export_37.txt"), sep = ",", header = TRUE)
                colnames(MART)[2] <- "Gene"
                
                normalizzati <- normalizzati %>%
                        left_join(MART, by = c("ID" = "Gene")) %>% 
                        dplyr::select(-ID) %>%
                        relocate(ID = Gene.stable.ID,  everything()) 
        }
        
        # Sort by variance and apply limit
        limit <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA[1])
        normalizzati_MOFA <- normalizzati %>%
                arrange(desc(variance)) %>%
                head(limit)
        
        # Save output file
        write.table(normalizzati_MOFA,
                    file = paste0(Cell_type, "_", args[1], "_vs_", args[2], "_normalized_vst_variance.tsv"),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        #================ PLOTTING VARIANCE ================
        pdf(file = paste0("Gene_expression_Variance_", Cell_type, "_", args[1], "_vs_", args[2], ".pdf"))
        
        ggplot(normalizzati, aes(x = variance)) +
                geom_histogram(fill = "lightblue", color = "darkblue", alpha = 0.7, bins = 1000) +
                stat_bin(bins = 1000) %>% print()
        
        ggplot(normalizzati, aes(x = variance)) +
                geom_density(color = "darkblue", fill = "lightblue") +
                geom_vline(aes(xintercept = 2.5)) %>% print()
        
        dev.off()
}



#STARTING STATISTICS ANALYSIS

# Load and prepare the transcriptomics data
load_and_prepare_data <- function(counts, Meta) {
        # Load data using vroom for faster data reading
        DGE3 <- vroom(counts, delim = "\t", col_names = TRUE)
        Metadata <- vroom(Meta, delim = "\t", col_names = TRUE)
        
        # Set row names of the DGE matrix to gene names
        rownames(DGE3) <- DGE3$ID
        
        # Remove the last column (gene names)
        #DGE3 <- DGE3[, -ncol(DGE3)]
        
        list(DGE3 = DGE3, Metadata = Metadata)
}


# Reorder DGE3 matrix columns based on Metadata$ID
reorder_data <- function(DGE3, Metadata) {
        tmp <- sapply(Metadata$ID, function(id) match(id, colnames(DGE3)))
        genes <- DGE3$ID
        DGE3 <- DGE3[, tmp]
        colnames(DGE3) <- Metadata$ID
        DGE3 <-as.data.frame(DGE3)
        rownames(DGE3) <- genes
        
        list(DGE3 = DGE3, Metadata = Metadata)
}


# Filter out 'neg' condition from DGE3 and Metadata
filter_condition <- function(DGE3, Metadata, comparison) {
        if(comparison == "neg"){
                out <- "pos"
        }else{out <- "neg"}
        num <- which(Metadata$condition == out)
        DGE3 <- DGE3 %>% dplyr::select(-all_of(num))
        Metadata <- Metadata[-num, ]
        
        list(DGE3 = DGE3, Metadata = Metadata)
}

# Function to prepare design formula
prepare_design_formula <- function(Metadata, Methods) {
        confounder_cols <- grep("_CONFOUNDER_", colnames(Metadata), value = TRUE)
        
        # Ensure correct column types (numeric or factor)
        for (col in confounder_cols) {
                if (grepl("_N$", col) && !is.numeric(Metadata[[col]])) {
                        Metadata[[col]] <- as.numeric(Metadata[[col]])
                }
                if (grepl("_F$", col) && !is.factor(Metadata[[col]])) {
                        Metadata[[col]] <- as.factor(Metadata[[col]])
                }
        }
        
        # Create the design formula dynamically
        covariates <- paste(confounder_cols, collapse = " + ")
        if (NORMALIZATION == "NO") {
                design_formula <- as.formula(paste("~", covariates, "+ CONDITION"))
        } else if (NORMALIZATION == "YES") {
                if (covariates == ""){
                design_formula <- as.formula(paste("~ 0 +", " CONDITION"))
                }else{ design_formula <- as.formula(paste("~ 0 +", covariates, "+ CONDITION"))}
        }
        
        list(design_formula = design_formula)
}


# Run DESeq2 analysis
run_deseq2_analysis <- function(DGE3, Metadata, Confounders_design, comparison) {
        library(DESeq2)
        
        dds <- DESeqDataSetFromMatrix(countData = DGE3,
                                      colData = Metadata,
                                      design = Confounders_design$design_formula)
        
        dds$condition <- relevel(as.factor(Metadata$CONDITION), ref = unlist(args[2]))
        dds <- DESeq(dds)
        
        res <- results(dds)
        
        resOrdered <- res[order(res$padj),]
        
        write.table(as.data.frame(resOrdered),
                    file= paste(args[1], comparison, "_", args[2],"_",Cell_type,".tsv" ,sep =""),sep = "\t",row.names = TRUE)
        
        Normalized_heatmap <- as.data.frame(vst(assay(dds))) #Normalization VST
        
        return(list(dds = dds, Normalized_heatmap = Normalized_heatmap))
}


# Run LIMMA analysis
run_limma_analysis <- function(DGE3, Metadata, Confounders_design, comparison) {
        library(limma)
        library(edgeR)
        
        # Set up design matrix
        design <- model.matrix(Confounders_design$design_formula, data = Metadata)
        
        contr.matrix <- makeContrasts(
                CTRLvsCONDITION = paste(colnames(design)[1], "-", colnames(design)[2]),
                levels = design
        )
        
        vfit <- lmFit(DGE3, design)
        vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
        efit <- eBayes(vfit, trend = TRUE)
        
        dds <- as.data.frame(topTreat(efit, coef = 1, n = Inf, adjust = "BH"))
        
        dds<-dds[,c(2,1,6,3,4,5)]
        colnames(dds) <- c("AveExpr", "log2FoldChange", "B" , "t", "pvalue", "padj"  )
        
        write.table(dds,
                    file= paste(args[1], comparison, "_", args[2],"_",Cell_type, ".tsv" ,sep =""),sep= "\t")
        
        Normalized_heatmap <- DGE3
        
        return(list(dds = dds, Normalized_heatmap = Normalized_heatmap))
}


#Common path for DESEQ2 and LIMMA


save_counts_and_metadata <- function(DGE3, Metadata, directory2, args, Cell_type) {
        dir.create(file.path(directory2, "Subpopulation_input"), showWarnings = TRUE, recursive = TRUE, mode = "0777")
        
        DGE3_modified <- DGE3
        DGE3_modified$ID <- rownames(DGE3_modified)
        DGE3_modified <- DGE3_modified %>% dplyr::select(ID, everything())
        
        write.table(DGE3_modified, 
                    file = file.path(directory2, "Subpopulation_input", paste0(args[1], "_pos_and_", args[2], "_", Cell_type, "_counts.tsv")),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        
        write.table(Metadata, 
                    file = file.path(directory2, "Subpopulation_input", paste0(args[1], "_pos_and_", args[2], "_", Cell_type, "_metadata.tsv")),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        
        gc()
}


annotate_genes <- function(DGE, Mart, GENE_ANNOTATION, directory2, args, Cell_type, comparison) {
        setwd(directory2)
        DGE <- read.table(file = paste(args[1], comparison, "_", args[2],"_", Cell_type,".tsv", sep =""), sep="\t", header = TRUE)
        DGE$X <- rownames(DGE)
        DGE <- DGE %>% dplyr::select(X, everything())
        
        if (GENE_ANNOTATION == "GENE_NAME") {
                DGE$Gene.name <- DGE$X
        } else {
                DGE <- merge(DGE, Mart, by.x = "X", by.y = "Gene.stable.ID")
        }
        return(DGE)
}


write_gct_file <- function(DGE3, dds, Mart, GENE_ANNOTATION, NORMALIZATION, args, directory2, Cell_type, comparison) {
        if (NORMALIZATION == "NO") {
                tab <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
        } else {
                tab <- DGE3
        }
        
        tab$X <- rownames(tab)
        
        if (GENE_ANNOTATION == "GENE_NAME") {
                tab <- tab %>% dplyr::select(X, everything())
        } else {
                tab <- merge(tab, Mart, by.x = "X", by.y = "Gene.stable.ID")
                tab <- tab[, -ncol(tab)]
        }
        
        colnames(tab)[1] <- "NAME"
        tab$Description <- "NA"
        tab <- tab[, c(1, ncol(tab), 3:(ncol(tab)-1))]
        
        gct_file <- file.path(directory2, paste0(args[1], comparison, "_", args[2], "_", Cell_type, "_RAW.gct"))
        if (file.exists(gct_file)) file.remove(gct_file)
        
        write("#1.2", file = gct_file)
        write(c(nrow(tab), ncol(tab) - 2), file = gct_file, append = TRUE, sep = "\t")
        write.table(tab, file = gct_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
}



save_metadata_files <- function(Metadata, args, Cell_type, directory2, comparison) {
        write.table(Metadata[order(Metadata$CONDITION, decreasing = TRUE), ],
                    file = file.path(directory2, paste0(args[1], comparison, "_", args[2], "_", Cell_type, "_RAW_meta.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
        write.table(Metadata[order(Metadata$CONDITION), c("ID", "CONDITION")],
                    file = file.path(directory2, paste0(args[1], comparison, "_", args[2], "_", Cell_type, "_RAW_meta_groups.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)
}


create_volcano_plot <- function(DGE, padju, log2FC, args, Cell_type, directory2, comparison) {
        
        library(ggplot2)
        library(ggrepel)
        
        ALTI <- subset(DGE, padj < padju & log2FoldChange > log2FC)
        ALTI <- ALTI[order(ALTI$padj), ]
        
        BASSI <- subset(DGE, padj < padju & log2FoldChange < -log2FC)
        BASSI <- BASSI[order(BASSI$padj), ]
        
        NO <- subset(DGE, padj >= padju)
        
        pdf(file = file.path(directory2, paste0("plot_DGE_", args[1], comparison, "_", args[2], "_", Cell_type, ".pdf")))
        
        p <- ggplot() +
                ggtitle(paste("DGE", Cell_type, args[1], comparison, "_", args[2])) +
                theme(plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5),
                      axis.title = element_text(size = 14, face = "bold")) +
                geom_point(data = ALTI, aes(x = log2FoldChange, y = -log10(padj)), color = "red") +
                geom_text_repel(data = head(ALTI, 25), aes(x = log2FoldChange, y = -log10(padj), label = Gene.name), size = 3) +
                geom_point(data = BASSI, aes(x = log2FoldChange, y = -log10(padj)), color = "blue") +
                geom_text_repel(data = head(BASSI, 25), aes(x = log2FoldChange, y = -log10(padj), label = Gene.name), size = 3) +
                geom_point(data = NO, aes(x = log2FoldChange, y = -log10(padj)), color = "black") +
                geom_vline(xintercept = c(-log2FC, log2FC), col = "red") +
                geom_hline(yintercept = -log10(padju), col = "red")
        
        print(p)
        dev.off()
        
        return(list(ALTI = ALTI, BASSI = BASSI))
}


write_gene_lists <- function(ALTI, BASSI, args, Cell_type, directory2, comparison) {
        write.table(ALTI$Gene.name, 
                    file = file.path(directory2, paste0(Cell_type, "_", args[1], comparison, "_", args[2], "_UP_ENRICHR.tsv")),
                    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        write.table(BASSI$Gene.name, 
                    file = file.path(directory2, paste0(Cell_type, "_", args[1], comparison, "_", args[2], "_DOWN_ENRICHR.tsv")),
                    sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        write.table(ALTI, 
                    file = file.path(directory2, paste0(Cell_type, "_", args[1], comparison, "_", args[2], "_GENESUP.tsv")),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        
        write.table(BASSI, 
                    file = file.path(directory2, paste0(Cell_type, "_", args[1], comparison, "_", args[2], "_GENESDOWN.tsv")),
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}


create_heatmap <- function(ALTI, BASSI, Normalized_heatmap, n_genes_heat, Mart, GENE_ANNOTATION, Metadata, Cell_type, args, directory2, COMMAND_ADVANCED) {
        if ((nrow(ALTI) + nrow(BASSI)) < 3) return(NULL)
        
        ALTI <- head(ALTI, n_genes_heat)
        BASSI <- head(BASSI, n_genes_heat)
        gene_ids <- c(ALTI$X, BASSI$X)
        Heat <- Normalized_heatmap[rownames(Normalized_heatmap) %in% gene_ids, ]
        rownames(Heat) <- gene_ids
        
        if (GENE_ANNOTATION == "ENSEMBL") {
                Heat$GENE <- rownames(Heat)
                Heat <- merge(Heat, Mart, by.x = "GENE", by.y = "Gene.stable.ID")
                Heat <- Heat[!duplicated(Heat[, ncol(Heat)]), ]
                rownames(Heat) <- Heat[, ncol(Heat)]
                Heat <- Heat[, -c(1, ncol(Heat))]
        }
        
        col_fun <- colorRamp2(c(min(Heat) - 3, mean(colMeans(Heat)), max(Heat) - 3), c("blue", "black", "yellow"))
        
        if("condition"%in% colnames(Metadata)){
        t <- setNames(c("blue", "red"), c(args[2], args[1]))
        t2 <- setNames(c("blue", "red", "green"), c("healthy", "pos", "neg"))
        ha <- HeatmapAnnotation(condition = Metadata$CONDITION, panel_status = Metadata$condition, col = list(condition = t, panel_status=t2))
        
        } else{
                
                t <- setNames(c("blue", "red"), c(args[2], args[1]))
                ha <- HeatmapAnnotation(condition = Metadata$CONDITION, col = list(condition = t))
                
        }
        
        pdf(file = file.path(directory2, paste0("Heatmap_top_genes_", Cell_type, "_", args[2], "_vs_", args[1], "_pos.pdf")), width = 10, height = 10)
        
        heat_obj <- Heatmap(
                Heat, km = 2, name = "vst_counts", col = col_fun,
                clustering_distance_rows = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[1],
                clustering_method_rows = "complete",
                clustering_method_columns = COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[2],
                row_dend_width = unit(0.5, "cm"),
                column_dend_height = unit(60, "mm"),
                column_names_gp = grid::gpar(fontsize = 6),
                row_names_gp = grid::gpar(fontsize = 8),
                top_annotation = ha
        )
        
        print(heat_obj)
        dev.off()
}


run_enrichr_analysis <- function(ALTI, BASSI, directory, args, Cell_type,directory_path, comparison) {
        
        library(enrichR)
        
        pdf(file = file.path(directory_path, "Pathways_DGE_EnrichR.pdf"), width = 20, height = 9)
        dbs <- c("Reactome_2022", "GO_Biological_Process_2023", "CODE_and_ChEA_Consensus_TFs_from_ChIP-X")
        
        if (getOption("enrichR.live")) {
                if (length(ALTI$Gene.name) > 0) {
                        enriched_up <- enrichr(ALTI$Gene.name, dbs)
                        for (i in seq_along(dbs)) {
                                if (nrow(enriched_up[[i]]) > 0) {
                                        write.table(enriched_up[[i]],
                                                    file = file.path(directory_path, "TABLES", paste0(names(enriched_up)[i], "_upregulated.tsv")),
                                                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
                                        print(plotEnrich(enriched_up[[i]], showTerms = 25, numChar = 40, y = "Count", orderBy = "P.value",
                                                         title = paste("Enrichment_analysis_upregulated_genes_", dbs[i])))
                                }
                        }
                }
                
                if (length(BASSI$Gene.name) > 0) {
                        enriched_down <- enrichr(BASSI$Gene.name, dbs)
                        for (i in seq_along(dbs)) {
                                if (nrow(enriched_down[[i]]) > 0) {
                                        write.table(enriched_down[[i]],
                                                    file = file.path(directory_path, "TABLES", paste0(names(enriched_down)[i], "_downregulated.tsv")),
                                                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
                                        print(plotEnrich(enriched_down[[i]], showTerms = 25, numChar = 40, y = "Count", orderBy = "P.value",
                                                         title = paste("Enrichment_analysis_downregulated_genes_", dbs[i])))
                                }
                        }
                }
        }
        
        dev.off()
}



