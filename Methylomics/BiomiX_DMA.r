#### INPUT PROMPT MANUAL (DEBUGGING) ----
# 
# args = as.list(c("Wholeblood","SLE"))
# print(args)
# 
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"mutated"
# args[2] <-"unmutated"
# args[3] <-"C:/Users/crist/Desktop/BiomiX2.5"
# 
# directory <- args[[3]]
# iterations = 1
# i=1
# selection_samples = "NO"
# Cell_type = "METHY"
# STATISTICS = "YES"
# DIR_METADATA <- readLines("C:/Users/crist/Desktop/BiomiX2.5/directory.txt")
# renv::load(paste(directory,"_INSTALL",sep="/"))

library(vroom)
library(dplyr)
library(readxl)
library(data.table)
library(matrixStats)
library(tidyverse)
library(rjson)



#### LOADING JSON FILE ----
setwd(directory)
combined_json <- jsonlite::fromJSON(txt = "COMBINED_COMMANDS.json")

directory2 <- paste(directory,"/Methylomics/",sep="")
setwd(directory2)

source(paste(directory2,"BiomiX_DMA_functions.r",sep=""))

directory2 <- paste(directory,"/Methylomics/INPUT",sep="")
setwd(directory2)


#External arguments required for the BiomiX_DMA.r script (methylomics)
COMMAND <- combined_json[["COMMANDS"]]
COMMAND_ADVANCED <- combined_json[["COMMANDS_ADVANCED"]]
LogFC <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METHYLOMICS[1])
padju <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METHYLOMICS[2])
array <- COMMAND_ADVANCED$ADVANCED_OPTION_METHYLOMICS[3]
n_genes_heat<-as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[3])
DIR_METADATA <- combined_json[["DIRECTORY_INFO"]][["METADATA_DIR"]]
Cell_type <- COMMAND$LABEL[i] #Cell_type is the label used to name saved files
args <- combined_json[["COMMAND_LINE_ARGS"]] #Belongs to the original arguments run in the BiomiX_beta.r
# i = iteration number from the BiomiX_beta.r Script
# STATISTICS = "YES" # Defined in the BiomiX_beta.r (If run the statistics on the omics data)




if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
        Metadata_total <- read_excel(DIR_METADATA)
        print("Metadata Excel File read successfully!")
}else{
        Metadata_total <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)}


if (grepl("\\.xlsx$|\\.xls$", COMMAND$DIRECTORIES[i])) {
        Matrix <- read_excel(COMMAND$DIRECTORIES[i])
        print("Metadata Excel File read successfully!")
}else{
        Matrix <- vroom(COMMAND$DIRECTORIES[i] , delim = "\t", col_names = TRUE)}



Metadata_Bcell <- Metadata_total %>%
  arrange(ID)

Identifier <- Matrix$ID
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

Metadata <- Metadata_Bcell




#FILTERING BASED ON METADATA

filter_metadata_matrix(COMMAND_ADVANCED, Metadata, Matrix) 


#Matrix and Metadata format rearrangement  
Metadata_individual <- Metadata

num<-Metadata_individual$ID[Metadata_individual$CONDITION == args[1]]
#Disease choice in metadata
numero <- Metadata_individual$ID[Metadata_individual$CONDITION == args[2]]
#Control choice in metadata
Matrix <- as.data.frame(Matrix)
Matrix$ID <- Identifier

Matrix <- Matrix[, c(ncol(Matrix), 1:(ncol(Matrix)-1))]

vector<- c(num,numero)
vector2 <- c("ID",vector)

Matrix<-Matrix[,colnames(Matrix) %in% vector2]
colnames(Matrix)[-1]<-colnames(Matrix)[-1][order(colnames(Matrix)[-1])]

Metadata_individual<-Metadata_individual[Metadata_individual$ID %in% vector,]
Metadata_individual<-Metadata_individual[order(Metadata_individual$ID),]

dim(Matrix)


if(COMMAND$PREVIEW[i] == "YES"){
        
#QC preview visualization

#Block to open the Quality control windows to select filtering, normalization and visualize the QC.
#ADD the ID to the first column
Samples_preview<-colnames(Matrix[,-1])
cpg <- Matrix$ID

numeric_data <- safe_transpose_df(Matrix[,-1])

rownames(numeric_data) <- Samples_preview
colnames(numeric_data) <- cpg

metadata <- Metadata_individual
numeric_data <- as.data.frame(numeric_data)
colnames(metadata)[colnames(metadata) == "CONDITION"] <- "condition"

# Source the BiomiX_preview script to load the runShinyApp() function
source(paste(directory,'/BiomiX_preview.r', sep=""))

browser_analysis <- readLines(paste(directory,'/_INSTALL/CHOISE_BROWSER_pre-view',sep=""), n = 1)

# Now call the runShinyApp function with numeric_data and metadata
options(browser = get_default_browser())
print("Pre QC data visualization. If the browser does not open automatically copy and paste the link on a browser")
print("One completed the analysis modification click on close app to continue the analysis")
print("ATTENTION!: IF the browser does not open set up the path to your browser on the CHOISE_BROWSER_pre-view file in the BiomiX _INSTALL folder")

Preview <- shiny::runApp(runShinyApp(numeric_data, metadata), launch.browser = TRUE)

Matrix<-as.data.frame(t(Preview$matrix))
Metadata_individual<-Preview$metadata
colnames(Metadata_individual)[colnames(Metadata_individual) == "condition"] <- "CONDITION"

Matrix <- cbind(cpg,Matrix)
colnames(Matrix)[1] <- "ID"


}else{
        print("no QC pre-visualization")
}


#### METHYLATION FILTER ON VARIANCE + INTEGRATION ----


varianze = NULL

Vari <- as.matrix(Matrix[,-1])

Vari <- as.data.frame(lapply(Vari, as.numeric))
varianze <- rowVars(as.matrix(Vari), na.rm = TRUE)

# varianze <-apply(Vari,1,var)
# for (i in 1:nrow(Vari)){
#         varianze <- append(varianze, var(as.numeric(Vari[i,])))
# }
gc()

Matrix2 <- Matrix
Matrix2$variance <- varianze

dir.create(path = paste(directory,"/Integration/INPUT/", "Methylome_",Cell_type,"_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory,"/Integration/INPUT/", "Methylome_",Cell_type, "_",args[1],"_vs_", args[2], sep ="") 

setwd(directory2)

#SELECT JUST THE TOP 10000 FEATURES WITH HIGHER VIABILITY
Matrix2<- Matrix2 %>% arrange(desc(variance))

limit<-as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA[1])
if(length(rownames(Matrix2)) > limit){
        Matrix2 <- Matrix2[1:limit,]}

write.table(x=Matrix2  , file= paste(Cell_type,"_matrix_MOFA.tsv",sep = "")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE)
write.table(x=Metadata_individual  , file= paste(Cell_type, "_metadata_MOFA.tsv", sep = "")  ,sep= "\t", row.names = F, col.names = T,  quote = FALSE,)


if (STATISTICS == "YES"){



#### STATISTICAL DMP ANALYSIS ----
  
  
  perform_statistical_dmp_analysis(
    Matrix = Matrix,
    Metadata_individual = Metadata_individual,
    directory = directory,
    Cell_type = Cell_type,
    args = args,
    array = array,
    padju = 0.05,
    LogFC = 0.25
  )
  
  
#HEATMAP

  
  create_heatmap(
    Matrix = Matrix,
    results = results,
    ALTI = ALTI,
    BASSI = BASSI,
    n_genes_heat = 50,
    Metadata_individual = Metadata_individual,
    directory2 = directory2,
    directory_path = directory_path,
    Cell_type = Cell_type,
    args = args
  )
  



#PATHWAY ANALYSIS
  perform_pathway_analysis(
    ALTI = ALTI,
    BASSI = BASSI,
    directory_path = directory_path
  )
  