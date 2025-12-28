
cat("\n\n\n\n\n          /////////      ///    /////////    ///       ///     ///    ///      ///\n         ///    ///     ///    ///   ///    /////  //////     ///      ///  ///  \n        /////////      ///    ///   ///    ///  ///  ///     ///        ///    \n       ///    ///     ///    ///   ///    ///  ///  ///     ///      ///  ///  \n      ///////////    ///    /////////    ///  ///  ///     ///    ///      ///\n                                                               ///            ///\n                                                                                 /////\n\n\n\n\n\n")

#directory <- "/media/henry/My_Passport/Article_MULTI_BETA"

# #MANUAL INPUT
# args = as.list(c("BLymphocytes","SLE"))
# args[1] <-"C1"
# args[2] <-"C2"
# args[3] <-"C:/Users/crist/Desktop/BiomiX2.5"



#Sys.sleep(5)
print("Welcome to BiomiX toolkit")

#Sys.sleep(5)

print("Information acquired, running the analysis")

#Sys.sleep(5)

#if(args[1] == ""){
        # taking input with showing the message
        args = commandArgs(trailingOnly=TRUE)
#}

#directory <- "/home/cristia/Scrivania/PhD/Bioinfo/Article_MULTI_BETA"
directory <- unlist(args[[3]])
setwd(paste(directory,"_INSTALL",sep="/"))
renv::load(paste(directory,"_INSTALL",sep="/"))

setwd(directory)
print(directory)

library(vroom)


#LOADING ARGUMENTS AND VARIABLES REQUIRED IN THE CODE

library(jsonlite)
library(tidyverse)
combined_json <- jsonlite::fromJSON(txt = "COMBINED_COMMANDS.json")


COMMAND <- combined_json[["COMMANDS"]]
COMMAND_MOFA <- combined_json[["COMMANDS_MOFA"]]
COMMAND_ADVANCED <- combined_json[["COMMANDS_ADVANCED"]]
DIR_METADATA <- combined_json[["DIRECTORY_INFO"]][["METADATA_DIR"]]
DIR_METADATA_output <- combined_json[["DIRECTORY_INFO"]][["OUTPUT_DIR"]]


#args = as.list(c("BLymphocytes","SJS"))
print(args)
#print(args[1])
#print(args[2])
#print(args[3])




#FIND A WAY TO SELECT ANY POSSIBLE INPUT WITHOUT CONSIDERING THE ORDER

n_iteration<- sum(COMMAND$DATA_TYPE == "Undefined")
position <- which(COMMAND$DATA_TYPE %in% "Undefined")
cat("\n\n\n\n\n Checking Undefined datasets.... \n\n\n\n\n")
iterations <- 0
for (i in position){
        if(COMMAND$ANALYSIS[i] == "YES" & COMMAND$INTEGRATION[i] == "YES"
           | COMMAND$ANALYSIS[i] == "YES" & COMMAND$INTEGRATION[i] == "NO")
           {
                
                STATISTICS <- "YES"
                Cell_type <- COMMAND$LABEL[i]
                cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
                iterations = iterations + 1
                selection_samples = COMMAND$SELECTION[i]
                #purity_filter =COMMAND$PURITY[i]
                directory2 <- paste(directory,"/Undefined/INPUT",sep="")
                #files <- grep("MATRIX*",list.files(directory2),value=TRUE)
                #files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
                source(paste(directory,"/Undefined/BiomiX_Undefined.r",sep=""))
                cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
                gc()
        }else if (COMMAND$ANALYSIS[i] == "NO" & COMMAND$INTEGRATION[i] == "YES") {
          
                #The following block allows you to re-run the single omics pipelines to generate the normalised data
                #Each time you run and integration analysis. It is quite time consuming, but gives stability.
          
                STATISTICS <- "NO"
                Cell_type <- COMMAND$LABEL[i]
                cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
                iterations = iterations + 1
                selection_samples = COMMAND$SELECTION[i]
                #purity_filter =COMMAND$PURITY[i]
                directory2 <- paste(directory,"/Undefined/INPUT",sep="")
                #files <- grep("MATRIX*",list.files(directory2),value=TRUE)
                #files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
                #source(paste(directory,"/Undefined/BiomiX_Undefined.r",sep=""))
                #cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
                cat("\n\n\n\n\n  Reminder: \n The Single Omics Analysis was not selected for", Cell_type, "\n If you are performing multiomics analysis and you did not run it before please do it,\n without the multiomics pipeline will not find the input data")
                
                gc()
                
        }
}




#===============================================================================






#FIND A WAY TO SELECT ANY POSSIBLE INPUT WITHOUT CONSIDERING THE ORDER

n_iteration<- sum(COMMAND$DATA_TYPE == "Transcriptomics")
position <- which(COMMAND$DATA_TYPE %in% "Transcriptomics")
cat("\n\n\n\n\nChecking Transcriptomics datasets.... \n\n\n\n\n")
iterations <- 0
for (i in position){
if(COMMAND$ANALYSIS[i] == "YES" & COMMAND$INTEGRATION[i] == "YES"
   | COMMAND$ANALYSIS[i] == "YES" & COMMAND$INTEGRATION[i] == "NO")
   {

        STATISTICS <- "YES"
        Cell_type <- COMMAND$LABEL[i]
        cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
        iterations = iterations + 1
        selection_samples = COMMAND$SELECTION[i]
        #purity_filter =COMMAND$PURITY[i]
        directory2 <- paste(directory,"/Transcriptomics/INPUT",sep="")
        #files <- grep("MATRIX*",list.files(directory2),value=TRUE)
        #files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
        source(paste(directory,"/Transcriptomics/Biomix_DGE_GENES_LIMMA.r",sep=""))
        cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
        gc()

        }else if (COMMAND$ANALYSIS[i] == "NO" & COMMAND$INTEGRATION[i] == "YES") {
        STATISTICS <- "NO"
        Cell_type <- COMMAND$LABEL[i]
        cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
        iterations = iterations + 1
        selection_samples = COMMAND$SELECTION[i]
        #purity_filter =COMMAND$PURITY[i]
        directory2 <- paste(directory,"/Transcriptomics/INPUT",sep="")
        #files <- grep("MATRIX*",list.files(directory2),value=TRUE)
        #files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
        #source(paste(directory,"/Transcriptomics/Biomix_DGE_GENES_LIMMA.r",sep=""))
        #cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
        cat("\n\n\n\n\n  Reminder: \n The Single Omics Analysis was not selected for", Cell_type, "\n If you are performing multiomics analysis and you did not run it before please do it,\n without the multiomics pipeline will not find the input data")
        gc()
}
}



#===============================================================================


#Sys.sleep(5)

n_iteration<- sum(COMMAND$DATA_TYPE == "Metabolomics")
position <- which(COMMAND$DATA_TYPE %in% "Metabolomics")
cat("\n\n\n\n\nChecking Metabolomics datasets.... \n\n\n\n\n")
iterations <- 0
for (i in position){
        if(COMMAND$ANALYSIS[i] == "YES" & COMMAND$INTEGRATION[i] == "YES"
           | COMMAND$ANALYSIS[i] == "YES" & COMMAND$INTEGRATION[i] == "NO")
           {
        
        STATISTICS <- "YES"
        Cell_type <- COMMAND$LABEL[i]
        cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
        iterations = iterations + 1
        selection_samples = COMMAND$SELECTION[i] # TO ADD
        #purity_filter =COMMAND$PURITY[i] #TO ADD
        directory2 <- paste(directory,"/Metabolomics/INPUT",sep="")
        files <- grep("MATRIX*",list.files(directory2),value=TRUE)
        files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
        source(paste(directory,"/Metabolomics/BiomiX_DMA.r",sep = ""))
        cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
        gc()
        
        }else if (COMMAND$ANALYSIS[i] == "NO" & COMMAND$INTEGRATION[i] == "YES") {
                STATISTICS <- "NO"
                Cell_type <- COMMAND$LABEL[i]
                cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
                iterations = iterations + 1
                selection_samples = COMMAND$SELECTION[i] # TO ADD
                #purity_filter =COMMAND$PURITY[i] #TO ADD
                directory2 <- paste(directory,"/Metabolomics/INPUT",sep="")
                files <- grep("MATRIX*",list.files(directory2),value=TRUE)
                files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
                #source(paste(directory,"/Metabolomics/BiomiX_DMA.r",sep = ""))
                #cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
                cat("\n\n\n\n\n  Reminder: \n The Single Omics Analysis was not selected for", Cell_type, "\n If you are performing multiomics analysis and you did not run it before please do it,\n without the multiomics pipeline will not find the input data")
                gc() 
        }
}

#================================================================================

n_iteration<- sum(COMMAND$DATA_TYPE == "Methylomics")
position <- which(COMMAND$DATA_TYPE %in% "Methylomics")
cat("\n\n\n\n\nChecking Methylomics datasets.... \n\n\n\n\n")
iterations <- 0
for (i in position){
        if(COMMAND$ANALYSIS[i] == "YES" & COMMAND$INTEGRATION[i] == "YES" 
           | COMMAND$ANALYSIS[i] == "YES" & COMMAND$INTEGRATION[i] == "NO")
                {
        
        STATISTICS <- "YES"
        Cell_type <- COMMAND$LABEL[i]
        cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
        iterations = iterations + 1
        selection_samples = COMMAND$SELECTION[i] # TO ADD
        #purity_filter =COMMAND$PURITY[i] #TO ADD
        directory2 <- paste(directory,"/Methylomics/INPUT",sep="")
        files <- grep("MATRIX*",list.files(directory2),value=TRUE)
        files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
        source(paste(directory,"/Methylomics/BiomiX_DMA.r",sep=""))
        cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
        gc()
        
        } else if (COMMAND$ANALYSIS[i] == "NO" & COMMAND$INTEGRATION[i] == "YES") {
                STATISTICS <- "NO"
                Cell_type <- COMMAND$LABEL[i]
                cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
                iterations = iterations + 1
                selection_samples = COMMAND$SELECTION[i] # TO ADD
                #purity_filter =COMMAND$PURITY[i] #TO ADD
                directory2 <- paste(directory,"/Methylomics/INPUT",sep="")
                files <- grep("MATRIX*",list.files(directory2),value=TRUE)
                files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
                #source(paste(directory,"/Methylomics/BiomiX_DMA.r",sep=""))
                #cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
                cat("\n\n\n\n\n  Reminder: \n The Single Omics Analysis was not selected for", Cell_type, "\n If you are performing multiomics analysis and you did not run it before please do it,\n without the multiomics pipeline will not find the input data")
                gc()
        }
}

#================================================================================


if(COMMAND_MOFA[2,2] == "YES") {

if(COMMAND_MOFA[1,2] == "MOFA_INTEGRATION") {        
        
        
        if(COMMAND_MOFA[3,2] == 0) {
                
                cat(paste( "\n\n\n\n\nStarting the AUTOMATIC MOFA analysis \n\n\n\n\n", sep =""))
                Sys.sleep(5)
                source(paste(directory,"/integration/MOFA_MULTI2_AUTO.R",sep=""))
                cat("\n\n\n\n\n AUTOMATIC MOFA analysis complete ^-^")
                
        }else{
                
                cat(paste( "\n\n\n\n\nStarting the MOFA analysis \n\n\n\n\n", sep =""))
                Sys.sleep(5)
                source(paste(directory,"/integration/MOFA_MULTI2.R",sep=""))
                cat("\n\n\n\n\n MOFA analysis complete ^-^")
        }
                
}
        if(COMMAND_MOFA[1,2] == "DIABLO_INTEGRATION") {
                
                cat(paste( "\n\n\n\n\nStarting the DIABLO analysis \n\n\n\n\n", sep =""))
                Sys.sleep(5)
                source(paste(directory,"/integration/DIABLO_MULTI2.R",sep=""))
                cat("\n\n\n\n\n DIABLO analysis complete ^-^")
                
        }
    
        if(COMMAND_MOFA[1,2] == "SNF_INTEGRATION") {
          
                cat(paste( "\n\n\n\n\nStarting the SNF analysis \n\n\n\n\n", sep =""))
                Sys.sleep(5)
                source(paste(directory,"/integration/SNF_int.R",sep=""))
                cat("\n\n\n\n\n SNF analysis complete ^-^")
          
        }
  
        if(COMMAND_MOFA[1,2] == "NEMO_INTEGRATION") {
    
                cat(paste( "\n\n\n\n\nStarting the NEMO analysis \n\n\n\n\n", sep =""))
                Sys.sleep(5)
                source(paste(directory,"/integration/NEMO_int.R",sep=""))
                cat("\n\n\n\n\n NEMO analysis complete ^-^")
    
  }
  
}


#================================================================================

cat(paste( "\n\n\n\n\nStarting the integrated factors correlation vs clinical factors \n\n\n\n\n", sep =""))

#Sys.sleep(5)

if(COMMAND_MOFA[1,2] == "DIABLO_INTEGRATION" | COMMAND_MOFA[1,2] == "MOFA_INTEGRATION" & COMMAND_MOFA[2,2] == "YES") {
        source(paste(directory,"/Clinical_data/BiomiX_Clinical.R",sep=""))
        cat("\n\n\n\n\n integrated factors correlation vs clinical factors completed ^-^")
}



#================================================================================

cat(paste( "\n\n\n\n\nStarting the integrated factors articles matching \n\n\n\n\n", sep =""))

#Sys.sleep(5)

if(COMMAND_MOFA[1,2] == "DIABLO_INTEGRATION" | COMMAND_MOFA[1,2] == "MOFA_INTEGRATION" & COMMAND_MOFA[2,2] == "YES") {
        source(paste(directory,"/integration/BiomiX_PUBMED.R",sep=""))
        cat("\n\n\n\n\n integrated factors articles matching completed ^-^")
}


#================================================================================

cat(paste( "\n\n\n\n\nStarting the integrated factors pathway analysis \n\n\n\n\n", sep =""))

#Sys.sleep(5)

if(COMMAND_MOFA[1,2] == "DIABLO_INTEGRATION" | COMMAND_MOFA[1,2] == "MOFA_INTEGRATION" & COMMAND_MOFA[2,2] == "YES") {
        source(paste(directory,"/integration/BiomiX_MOFA_MINER.R",sep=""))
        cat("\n\n\n\n\n integrated factors pathway analysis completed ^-^")
}


#================================================================================
#args[1] = "ICPLUS"

cat(paste( "\n\n\n\n\n Saving the files in the directory \n\n\n\n\n", sep =""))

Sys.sleep(5)

#directory2 <- paste(directory, "/",COMMAND$DATA_TYPE[1], "/OUTPUT/",COMMAND$LABEL[1], "_",args[1], "_vs_", args[2], sep ="")
#print(directory2)
#setwd(directory2)

#COPY OMICS DIRECTORY

if (directory != args[3]){

for(output in 1:length(na.omit(COMMAND$LABEL))) {
        print(output)
        print(paste(directory,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))
        if(file.exists(paste(directory,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))){
                dir.create(path =  paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                file.copy(paste(directory,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""), paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
                print(paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))
}
}

if(COMMAND_MOFA[2,2] == "YES") {

#COPY MOFA FILES
output<-list.files(path = paste(directory,"/","MOFA", "/OUTPUT/",sep =""), pattern = paste("^","AUTOMATIC_SEARCH_BEST_FACTORS_",args[1], "_vs_", args[2], sep=""))
if(file.exists(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""))){
        dir.create(path =  paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        file.copy(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep =""), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE) 
        print(paste(output,"copied", sep = " "))
}

#COPY MOFA DIRECTORY
file_to_transfer<-dir(path = paste(directory,"/","MOFA", "/OUTPUT/",sep =""), pattern = paste("^","MOFA_",args[1], "_vs_", args[2], sep="") )
for(output in file_to_transfer){
        if(file.exists(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""))){
                dir.create(path =  paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                file.copy(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
                #unlink(paste(directory,"/","MOFA", "/OUTPUT/",output, sep =""), recursive=TRUE)
                print(paste(output,"copied", sep = " "))
        }}

if(COMMAND_ADVANCED[1,12] == "YES" | COMMAND_ADVANCED[2,12]  == "YES") {

output<-dir(path = paste(directory,"/","Clinical_data", "/OUTPUT/",sep =""), pattern = paste("^","Clinical_",args[1], "_vs_", args[2], sep="") )
        if(file.exists(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""))){
                dir.create(path =  paste(DIR_METADATA_output,"/","Clinical_data","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                file.copy(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","Clinical_data","/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
                #unlink(paste(directory,"/","Clinical_data", "/OUTPUT/",output, sep =""), recursive=TRUE)
                print(paste(output,"copied", sep = " "))
        }
}
}
}

cat("\n\n\n\n\n file saving complete ^-^")



# ---------- HELPERS ----------

make_param_table <- function(names, values) {
  data.frame(
    Parameter = names,
    Value = values,
    stringsAsFactors = FALSE
  )
}

write_section <- function(section_title, table = NULL) {
  write(section_title, file = out_file, append = TRUE)
  if (!is.null(table)) {
    write.table(
      table, file = out_file, append = TRUE,
      row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE
    )
  }
}

# ---------- BUILD PARAMETER TABLES ----------

matr <- make_param_table(
  names = c(
    "LOG2FC_TRANSCRIPTOMICS", "P.ADJ_TRANSCRIPTOMICS", "GENE_PANEL_FILE",
    "LOG2FC_METABOLOMICS", "P.ADJ_METABOLOMICS", "CPU_USED",
    "LOG2FC_METHYLOMICS", "P.ADJ_METHYLOMICS", "ARRAY_TYPE",
    "PANEL POSITIVITY WITH N°GENES WITHIN THE PANEL WITH Z-SCORE > 2",
    "PANEL POSITIVITY WITH N°GENES WITHIN THE PANEL WITH Z-SCORE > 1",
    "REMOVE CONTROL POSITIVE FOR GENE PANEL?",
    "CLUSTERING DISTANCE", "CLUSTERING METHOD", "N° GENES TO VISUALIZE IN THE HEATMAP?",
    "N° MOFA INPUT FEATURES"
  ),
  values = c(
    COMMAND_ADVANCED$ADVANCED_OPTION_TRASCRIPTOMICS,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS,
    COMMAND_ADVANCED$ADVANCED_OPTION_METHYLOMICS,
    COMMAND_ADVANCED$ADVANCED_OPTION_SUBGROUPING,
    COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS,
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[3]
  )
)

matr_2 <- make_param_table(
  names = c(
    "LEVEL_ANNOTATION_METABOLOMICS_DATASET?",
    "IF_ANNOTATED_WHICH_TYPE_OF_ANNOTATION?",
    "ION_MODE_ONLY_MS1_MODE", "NEURAL_MODE_ONLY_MS1_MODE", "TOLERANCE_IN_PPM_ONLY_MS1_MODE",
    "POSITIVE_ADDUCTS_ONLY_MS1_MODE", "NEGATIVE_ADDUCTS_ONLY_MS1_MODE", "DATABASE_FOR_MS1_ANNOTATION_ONLY_MS1_MODE",
    "INDEXING_1", "INDEXING_2", "INDEXING_3",
    "FILE_MS1_ANNOTATION_INDEXING_1_ONLY_MS1_MODE",
    "FILE_MS1_ANNOTATION_INDEXING_2_ONLY_MS1_MODE",
    "FILE_MS1_ANNOTATION_INDEXING_3_ONLY_MS1_MODE",
    "TOLERANCE_PPM_MATCH_MS1_MS2_PEAKS", "RETENTION_TIME_MATCH_MS1_MS2",
    "MS2_DATABASE_ANNOTATION_PRIORITY",
    "ION_MODE_ONLY_MS1-2_MODE",
    "POSITIVE_ADDUCTS__MS1-2_MODE", "NEGATIVE_ADDUCTS__MS1-2_MODE", "TYPE_LC_COLUMN",
    "DATABASE_FOR_MS1_ANNOTATION_ONLY_MS1-2_MODE", "TOLERANCE_IN_PPM_ONLY_MS1-2_MODE",
    "FILE_MS1_ANNOTATION_INDEXING_1_ONLY_MS1-2_MODE",
    "FILE_MS1_ANNOTATION_INDEXING_2_ONLY_MS1-2_MODE",
    "FILE_MS1_ANNOTATION_INDEXING_3_ONLY_MS1-2_MODE",
    "FILE_MS1_ANNOTATION_INDEXING_1_ONLY_MS1-2_MODE",
    "FILE_MS1_ANNOTATION_INDEXING_2_ONLY_MS1-2_MODE",
    "FILE_MS1_ANNOTATION_INDEXING_3_ONLY_MS1-2_MODE",
    "DIRECTORY_TO_MS2_FILES"
  ),
  values = c(
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3[1],
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY[1]
  )
)

matr_3 <- make_param_table(
  names = c(
    "FILTERING_METADATA_1", "TYPE_OF_DATA_1", "THRESHOLD_OR_NAME_1",
    "FILTERING_METADATA_2", "TYPE_OF_DATA_2", "THRESHOLD_OR_NAME_2",
    "FILTERING_METADATA_3", "TYPE_OF_DATA_3", "THRESHOLD_OR_NAME_3",
    "FILTERING_METADATA_4", "TYPE_OF_DATA_4", "THRESHOLD_OR_NAME_4"
  ),
  values = c(
    COMMAND_ADVANCED$ADVANCED_OPTION_METADATA_FILTERING_1,
    COMMAND_ADVANCED$ADVANCED_OPTION_METADATA_FILTERING_2,
    COMMAND_ADVANCED$ADVANCED_OPTION_METADATA_FILTERING_3,
    COMMAND_ADVANCED$ADVANCED_OPTION_METADATA_FILTERING_4
  )
)

matr_4 <- make_param_table(
  names = c(
    "MAX_ITERATION", "CONVERGENCE_SPEED", "THRESHOLD_CONTRIBUTION_WEIGHT",
    "TYPE_OF_RESEARCH", "N°_OF_ARTICLES",
    "N°_TOP_CONTRIBUTORS_TO_SEARCH_IN_PUBMED",
    "N°_KEYWORDS_EXTRACTED_BY_THE_KEYWORD_GENERATOR",
    "P.ADJ_THRESHOLD_PATHWAY", "N°_OF_PATHWAY_TO_VISUALIZE",
    "NUMERIC_CLINICAL_DATA_UPLOADED?", "BINARY_CLINICAL_DATA_UPLOADED)",
    "NUMERIC_CLINICAL_DATA_FILE?", "BINARY_CLINICAL_DATA_FILE",
    "SNF_N°Neighbours", "SNF_Variance_affinity_Matrix","SNF_N°Iteration",
    "NEMO_N°Neighbours", "NEMO_Variance_affinity_Matrix",
    "SNF_NEMO_Enrichment_Variable", "SNF_NEMO_Survival_Variable", "SNF_NEMO_Ground_truth_variable",
    "SNF_NEMO_Max_cluster_expected", "SNF_NEMO_N°feature_heatmap", "SNF_NEMO_N°input_features",
    "DIABLO_Design_matrix_path", "Diablo_features_selection?", "Diablo_N°iterations"
    
  ),
  values = c(
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA,
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY_2[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL,
    COMMAND_ADVANCED$ADVANCED_OPTION_CLINIC_DATA_DIRECTORY,
    COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS,
    COMMAND_ADVANCED$ADVANCED_OPTION_NEMO_OPTIONS,
    COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS,
    COMMAND_ADVANCED$ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN[1],
    COMMAND_ADVANCED$ADVANCED_OPTION_DIABLO_OPTIONS[1:2]
  )
)

# ---------- OUTPUT DIRECTORY ----------

dir_out <- file.path(directory, "Report_parameters")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)


library(openxlsx)

excel_file <- file.path(
  dir_out,
  paste0(
    "Report_", args[1], "_vs_", args[2], "_",
    format(Sys.time(), "%Y%m%d_%H%M"),
    ".xlsx"
  )
)

wb <- createWorkbook()

addWorksheet(wb, "OMICS_COMMANDS")
writeData(wb, "OMICS_COMMANDS", COMMAND)

addWorksheet(wb, "MOFA_COMMANDS")
writeData(wb, "MOFA_COMMANDS", COMMAND_MOFA)

addWorksheet(wb, "ADV_GENERAL")
writeData(wb, "ADV_GENERAL", matr)

addWorksheet(wb, "ADV_METABOLOMICS")
writeData(wb, "ADV_METABOLOMICS", matr_2)

addWorksheet(wb, "ADV_METADATA")
writeData(wb, "ADV_METADATA", matr_3)

addWorksheet(wb, "ADV_MOFA")
writeData(wb, "ADV_MOFA", matr_4)

saveWorkbook(wb, excel_file, overwrite = TRUE)
