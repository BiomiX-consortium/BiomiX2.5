#Remove selection from Biomix interface
#ADD variables for the correction of the linear model (no N_CONFOUNDER)




#### INPUT PROMPT ----
# 
# 
# MANUAL INPUT
# library(vroom)
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"mutated"
# args[2] <-"unmutated"
# args[3] <-"C:/Users/crist/Desktop/BiomiX2.5"
# #
# directory <- args[3]
# iterations = 1
# i=1
# selection_samples = "NO"
# purity_filter = "NO"
# Cell_type = "RNA"
# setwd(directory[[1]])
# COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
# COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")
# DIR_METADATA <- readLines("directory.txt")
# STATISTICS = "YES"
# renv::load(paste(directory,"_INSTALL",sep="/"))


#### DATASET REARRANGEMENT ----

print("starting data upload and preparation")

library(vroom)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(readxl)

directory2 <- paste(directory,"/Transcriptomics/",sep="")
setwd(directory2)

source(paste(directory2,"Biomix_DGE_GENES_LIMMA_FUNCTIONS.r",sep=""))

directory2 <- paste(directory,"/Transcriptomics/INPUT",sep="")
setwd(directory2)

print(iterations)

COMMAND$DIRECTORIES[i]
if (grepl("\\.xlsx$|\\.xls$", COMMAND$DIRECTORIES[i])) {
        Matrix <- read_excel(COMMAND$DIRECTORIES[i])
        print("Metadata Excel File read successfully!")
        }else{
                Matrix <-vroom(COMMAND$DIRECTORIES[i] , delim = "\t", col_names = TRUE)}


if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
        Metadata_total <- read_excel(DIR_METADATA)
        print("Metadata Excel File read successfully!")
        }else{
                Metadata_total <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)}

COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
log2FC <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[1])
padju <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[2])
Gene_panel <- COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3]
n_genes_heat<-as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[3])



#Debug
# Gene_panel <- "GENE_PANEL.txt"
# COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3] = "YES"
# #Conversion gene name
# Mart<-read.table("/home/cristia/BiomiX2.4/MOFA/x_BiomiX_DATABASE/mart_export_37.txt", sep = ",", header=TRUE)
# Matrix <- merge(Matrix, Mart, by.x = "ID", by.y = "Gene.stable.ID")
# Matrix$ID <- Matrix$Gene.name
# Matrix <- Matrix[,-ncol(Matrix)]
#==============


genes <- Matrix$ID



#TEST IF ALL MATRIX SAMPLES ARE INCLUDED IN METADATA

Metadata_Bcell <- Metadata_total %>%
        arrange(ID)

Matrix <- Matrix[, -1, drop = FALSE] %>% as.matrix()
Matrix <- Matrix[, order(colnames(Matrix)), drop = FALSE]

Metadata_Bcell <- Metadata_total
num <- which(colnames(Matrix) %in% Metadata_Bcell$ID)
Matrix <- Matrix[, num, drop = FALSE]

num <- which(Metadata_Bcell$ID %in% colnames(Matrix))
Metadata_Bcell <- Metadata_Bcell[num, ] %>%
        distinct(ID, .keep_all = TRUE) %>%
        arrange(ID)

tryCatch({
        Matrix <- Matrix[, order(colnames(Matrix)), drop = FALSE]
}, error = function(e) {
        message("Oops! An error occurred: ", e$message)
        message("Are you sure the transcriptomics matrix samples are the same as the metadata file?")
})


rownames(Matrix) <- genes
summary(as.factor(Metadata_Bcell$CONDITION))

Metadata_total = NULL




# Check type of gene name used (ENSEMBL or GENE SYMBOL)

Mart<-read.table(paste(directory,"/Integration/x_BiomiX_DATABASE/mart_export_37.txt", sep=""), sep = ",", header=TRUE)

result <- process_gene_annotation(Matrix, genes, directory)

DGE2 <- result$DGE2
genes_name <- result$genes_name
genes_name_C <- result$genes_name_C
GENE_ANNOTATION <- result$GENE_ANNOTATION

#Order Matrix samples and metadata
DGE2 <- DGE2[,order(colnames(DGE2))]
Metadata_Bcell <- arrange(Metadata_Bcell,ID)

#Chosen disease
num <- Metadata_Bcell$CONDITION == args[1]
#Choise control
numero <- Metadata_Bcell$CONDITION == args[2]

Metadata_Bcell <- Metadata_Bcell[num | numero,]
DGE2 <- DGE2[, num | numero]

if(any(num) & any(num)){ print("Conditions detected in Metadata samples")
        }else{print("ERROR: Conditions undetected in Metadata samples")}


#FILTERING SAMPLES BASED ON METADATA CRITERIA SELECTED. 

result <- filter_metadata_and_expression(Metadata_Bcell, DGE2, COMMAND_ADVANCED, directory)

Metadata_Bcell <- result$Metadata_Bcell
DGE2 <- result$DGE2

Metadata_Bcell$CONDITION
summary(as.factor(Metadata_Bcell$CONDITION))


#Preview visualization

if(COMMAND$PREVIEW[i] == "YES"){
        source(paste(directory,'/BiomiX_preview.r', sep=""))
        browser_analysis <- readLines(paste(directory,'/_INSTALL/CHOISE_BROWSER_pre-view',sep=""), n = 1)
        options(browser = get_default_browser())
        print(get_default_browser())
        result <- launch_qc_preview(DGE2, Metadata_Bcell, directory, browser_analysis)
        DGE2 <- result$DGE2
        Metadata_Bcell <- result$Metadata_Bcell
} else {
        print("no QC pre-visualization")
}

print("samples post pre-visualization ")
summary(as.factor(Metadata_Bcell$CONDITION))

#STATEMENT IF DATA ARE NORMALIZED OR NOT, IT LOOKS IF NUMERS ARE FLOATS OR INTEGER

results <- handle_normalization_and_export(DGE2, Metadata_Bcell, Cell_type, args, directory)

NORMALIZATION <- results$NORMALIZATION
DGE3 <- results$DGE3
DGE2 <- results$DGE2

if(!is.null(results$dds)){ 
        dds <- results$dds }else{ dds <- NULL}
   

#### HEATMAP SECTION ----

Metadata_Bcell <- generate_heatmap_signature(
        DGE2 = DGE2,
        dds = dds,
        Metadata_Bcell = Metadata_Bcell,
        Gene_panel = Gene_panel,
        GENE_ANNOTATION = GENE_ANNOTATION,
        Mart = Mart,
        directory = directory,
        Cell_type = Cell_type,
        args = args,
        COMMAND_ADVANCED = COMMAND_ADVANCED
)

  
  
  #### SAVING NORMALIZED DATA ----
  
  
  # Main refactored logic
  if (COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3] != "X") {
          Metadata <- prepare_metadata_conditions(Metadata_Bcell)
          #Filtering step, if yes the CTRL positive are removed otherwise the are not.
          if (COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING[3] != "YES") {
                  Metadata_all <- Metadata
                  DGE2 <- format_expression_matrix(DGE2, Metadata, genes_name, NORMALIZATION)
          } else {
                  filtered <- filter_samples_by_condition(DGE2, DGE3, Metadata)
                  DGE2 <- filtered$DGE2
                  DGE3 <- filtered$DGE3
                  Metadata <- filtered$Metadata
                  Metadata_all <- Metadata
                  DGE2 <- format_expression_matrix(DGE2, Metadata, genes_name, NORMALIZATION)
          }
          
  } else {
          Metadata <- Metadata_Bcell
          DGE2 <- as.data.frame(t(DGE2))
          # if (NORMALIZATION == "NO") {
          #         DGE2 <- DGE2[-nrow(DGE2), ]
          # }
          Metadata_all <- Metadata
          # Optional: uncomment if you want to keep formatting consistent
          # colnames(DGE2) <- c(genes_name)
          # rownames(DGE2) <- Metadata$ID
          # DGE2$samples <- Metadata$ID
  }
  
  # Save the outputs
  save_outputs(Metadata, DGE3, directory2, Cell_type, args[1], args[2])
  
  # Store output file names (optional use later)
  counts <- paste("gene_count_", Cell_type, "_", args[1], "_vs_", args[2], "_unnormalized.tsv", sep = "")
  Meta <- paste("Metadata_", Cell_type, "_", args[1], "_vs_", args[2], "_final.tsv", sep = "")
  
  
 
#### GENERATION MATRIX FOR DATA INTEGRATION ----
  
  generate_mofa_input(directory, Cell_type, args, NORMALIZATION, GENE_ANNOTATION, COMMAND_ADVANCED)
  

if (STATISTICS == "YES"){
        
        
#STARTING STATISTICS ANALYSIS
        
                dim(DGE3)
                dim(Metadata)
                # Load and prepare data
                setwd(directory2)
                data <- load_and_prepare_data(counts, Meta)
                DGE3 <- data$DGE3
                Metadata <- data$Metadata
                
                # Reorder data based on metadata
                data_backup <- reorder_data(DGE3, Metadata)
                DGE3 <- data_backup$DGE3
                Metadata <- data_backup$Metadata
                
                if (COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS[3] != "X") {
                        #Panel_state defines the ending part of the folders and files
                        #If the gene panel is set up it will include pos and neg comparison, otherwise only the CTRL vs DISEASE condition "" = EMPTY.
                        Panel_state <- c("", "pos", "neg") } else { Panel_state <- c("") }
                
                
for (comparison in Panel_state){
        
        DGE3 <- data_backup$DGE3
        Metadata <- data_backup$Metadata
                        
        if (comparison != ""){
                
                if (sum(Metadata$condition == comparison) > 2){
                
                # Filter out "neg" condition
                data <- filter_condition(DGE3, Metadata, comparison)
                DGE3 <- data$DGE3
                Metadata <- data$Metadata 
                
                #create directory DGE
                dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1], "_vs_", args[2],"/",args[1] ,comparison,"_vs_",args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                directory2 <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1], "_vs_", args[2],"/",args[1],comparison,"_vs_",args[2], sep ="")
                setwd(directory2)
                
                #create directory EnrichR
                directory_path <- file.path(directory, "Transcriptomics/OUTPUT", paste0(Cell_type, "_", args[1], "_vs_", args[2]),
                                            paste0(args[1], comparison, "_vs_", args[2]), "Pathway_analysis")
                dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
                dir.create(file.path(directory_path, "TABLES"), showWarnings = FALSE)
                
                } else{
                        print(paste(comparison,"samples lower than 3, the DGE analysis is not possible"))
                        next
                        }
                        
        } else {
                
                #create directory DGE
                dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                dir.create(path = paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1],"_vs_",args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                directory2 <- paste(directory,"/Transcriptomics/OUTPUT/", Cell_type, "_",args[1],"_vs_", args[2],"/",args[1], "_vs_",args[2], sep ="")
                setwd(directory2)
                
                #create directory EnrichR
                directory_path <- file.path(directory, "Transcriptomics/OUTPUT", paste0(Cell_type, "_", args[1], "_vs_", args[2]),
                                            paste0(args[1], "_vs_", args[2]), "Pathway_analysis")
                dir.create(directory_path, recursive = TRUE, showWarnings = FALSE)
                dir.create(file.path(directory_path, "TABLES"), showWarnings = FALSE)
                
                DGE3 <- DGE3
                Metadata <- Metadata 
                print("Control-disease comparison")
        }
                
        
                # Prepare design formula
                Confounders_design <- prepare_design_formula(Metadata, NORMALIZATION)  # or "LIMMA"
                
                
                # Run analysis based on normalization
                if (NORMALIZATION == "YES") {
                        dds <- run_limma_analysis(DGE3, Metadata, Confounders_design, comparison )
                } else {
                        dds <- run_deseq2_analysis(DGE3, Metadata, Confounders_design, comparison)
                }
                
                # Common part for DESEQ and LIMMA (MERGIN BRANCH)
                Normalized_heatmap <- dds$Normalized_heatmap
                dds <- dds$dds
                save_counts_and_metadata(DGE3, Metadata, directory2, args, Cell_type)
                DGE<-annotate_genes(DGE, Mart, GENE_ANNOTATION, directory2, args, Cell_type,comparison)
                write_gct_file(DGE3, dds, Mart, GENE_ANNOTATION, NORMALIZATION, args, directory2, Cell_type, comparison)
                save_metadata_files(Metadata, args, Cell_type, directory2, comparison)
                
                #Volcanoplot
                UP_DOWN_genes<- create_volcano_plot(DGE, padju, log2FC, args, Cell_type, directory2, comparison)
                #Collection upregulated genes in ALTI variable
                ALTI <- UP_DOWN_genes$ALTI
                #Collection downregulated genes in BASSI variable
                BASSI <- UP_DOWN_genes$BASSI
                write_gene_lists(ALTI, BASSI, args, Cell_type, directory2,comparison)
                create_heatmap(ALTI, BASSI, Normalized_heatmap, n_genes_heat, Mart, GENE_ANNOTATION, Metadata, Cell_type, args, directory2, COMMAND_ADVANCED)
                
                #EnrichR
                run_enrichr_analysis(ALTI, BASSI, directory, args, Cell_type, directory_path, comparison)
                
                }
                
                
                        }else{
        print("No statistical analysis")
}

gc()


