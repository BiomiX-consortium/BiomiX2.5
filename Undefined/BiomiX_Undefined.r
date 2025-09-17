library(vroom)
library(dplyr)
library(stringr)
library(rlist)
library(tibble)
library(readxl)

# MANUAL INPUT
# # # #
# library(vroom)
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"mutated"
# args[2] <-"unmutated"
# args[3] <-"C:/Users/crist/Desktop/BiomiX2.5"
# 
# directory <- args[[3]]
# selection_samples = "NO"
# Cell_type = "populations"
# DIR_METADATA <- readLines("C:/Users/crist/Desktop/BiomiX2.5/directory.txt")
# i = 1
# renv::load(paste(directory,"_INSTALL",sep="/"))



setwd(directory)

#LOADING ARGUMENTS AND VARIABLES REQUIRED IN THE CODE


library(jsonlite)
library(tidyverse)
combined_json <- jsonlite::fromJSON(txt = "COMBINED_COMMANDS.json")

COMMAND <- combined_json[["COMMANDS"]]
COMMAND_ADVANCED <- combined_json[["COMMANDS_ADVANCED"]]
COMMAND_MOFA <- combined_json[["COMMANDS_MOFA"]]
DIR_METADATA <- combined_json[["DIRECTORY_INFO"]][["METADATA_DIR"]]
DIR_METADATA_output <- combined_json[["DIRECTORY_INFO"]][["OUTPUT_DIR"]]

Heatmap_genes <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[3])

LogFC <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS[1])
padju <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS[2])


directory2 <- paste(directory,"/Undefined/",sep="")
setwd(directory2)

source(paste(directory2,"BiomiX_Undefined_functions.r",sep=""))

#loading matrix
directory2 <- paste(directory,"/Undefined/INPUT",sep="")
setwd(directory2)



print(i)
print(COMMAND$LABEL[i])

if (grepl("\\.xlsx$|\\.xls$", COMMAND$DIRECTORIES[i])) {
        matrix <- read_excel(COMMAND$DIRECTORIES[i])
        print("Matrix Excel File read successfully!")
}else{
        matrix <- vroom(COMMAND$DIRECTORIES[i],delim="\t",col_names = TRUE, comment = "#")
}


sam <- colnames(matrix)
pea <-matrix$ID
matrix <- t(matrix[,-1])
matrix <- add_column(as.data.frame(matrix), sam[-1], .after = 0) 
colnames(matrix) <-c("ID",pea)


#create directory
dir.create(path = paste(directory,"/Integration/INPUT/", "Undefined_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory, "/Integration/INPUT/", "Undefined_", Cell_type,  "_",args[1],"_vs_", args[2], sep ="")

setwd(directory2)


if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
        Metadata_total <- read_excel(DIR_METADATA)
        print("Metadata Excel File read successfully!")
}else{
        Metadata_total <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)}


#If statement in case of sample selection


        
process_metadata_matrix(matrix, Metadata_total)
        

#The metadata file contains the metadata of only the samples within the omics
#Excluding the ones not included in the original metadata file.
Metadata <- Metadata_Bcell
Metadata_individual=NULL
Metadata_reads=NULL
Metadata_Bcell=NULL
Metadata_tot = NULL
Metadata_total = NULL

#Check if metadata and matrix contain the same samples.
print("Same samples in matrix and metadata?")
print(all(Metadata$ID == matrix[,1]))
outs <- Metadata$CONDITION  == args[2]|Metadata$CONDITION  == args[1]
Metadata<- Metadata[outs,]
matrix<- matrix[outs,]


#FILTERING BASED ON METADATA
#Based on the criteria defined in the advance options.
apply_metadata_filter(Metadata, matrix, COMMAND_ADVANCED)



###
matrix<- as.data.frame(matrix)
rownames(matrix) <- matrix$ID
###


#This part of code is to refill the column having an NA name and rename the duplicated ones
na_cols <- which(is.na(colnames(matrix)))
new_names <- paste("UNKNOWN",1:sum(is.na(colnames(matrix))),sep="") 
colnames(matrix)[na_cols] <- new_names

colnames(matrix) <- make.unique(names(matrix), sep = "_")

#====================================
if(COMMAND$PREVIEW[i] == "YES"){
        
        #Block to open the Quality control windows to select filtering, normalization and visualize the QC.
        #ADD the ID to the first column
        browser_analysis <- readLines(paste(directory, '/_INSTALL/CHOISE_BROWSER_pre-view', sep = ""), n = 1)
        run_qc_preview(matrix, Metadata, directory)
        
}else{
        print("no QC pre-visualization")
}


#====================================
matrixs <- matrix 


matrixs <- add_column(matrixs, Metadata$CONDITION, .after = 1)
colnames(matrixs)[2] <- "CONDITION"


tryCatch({
        # If there is a missmatch between metadata and matrix this line will not work
        matrixs[,3:ncol(matrixs)]<-apply(matrixs[,3:ncol(matrixs)],2,as.numeric)
        
}, error = function(e) {
        # Personalized error message
        print(paste("Oops! An error occurred:", e$message))
        print("Are you sure the Undefined matrix samples are the same of the metadata file?")
})



#Calculus variance
varianze = NULL
Vari <- as.matrix(matrixs[,c(-1:-2)])
Vari <- as.data.frame(Vari)

#generazione e filtraggio per varianza
for (i in 1:ncol(Vari)){
        varianze <- append(varianze, var(as.numeric(Vari[,i])))
}


normalizzati<- as.data.frame(matrixs)
variance <- log10(varianze)

ordering<-data.frame(colnames(Vari), variance)
ordering$variance<-abs(ordering$variance)
ordering <- ordering[order(ordering$variance, decreasing = TRUE), ]
limit<-as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA[1])

#SAVING THE INPUT FOR THE MOFA ANALYSIS
#These data should not be normalized by log.

if(length(rownames(ordering)) > limit){
        ordering<-ordering[!ordering$variance == Inf,]
        normalizzati_MOFA <- ordering[1:limit,]
        ond<-colnames(matrixs) %in% normalizzati_MOFA[,1] 
        ond[1:2] <- TRUE
        normalizzati_MOFA<-matrixs[,ond]
}else{
normalizzati_MOFA <- matrixs}

write.table(normalizzati_MOFA,paste(directory2,"/Undefined_",Cell_type, "_MOFA.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")


if (STATISTICS == "YES"){
        




# #### DMS DIFFERENTIAL METABOLITE SIGNALS ANALYSIS----

CON <-DMS_REFERENCE(args[2])        
TEST <- DMS_TESTED(args[1])        
TEST$NAME <- row.names(TEST)
gc()

#Saving results
        
run_differential_analysis(TEST, matrix, Metadata, directory, Cell_type, args, padju, LogFC)
        
        
        ###  HEATMAP ###
        
run_differential_analysis(TEST, matrix, Metadata, directory, Cell_type, args, padju, LogFC)

        
}else{
        print("No statistical analysis and annotation")
}                
        


