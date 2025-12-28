


# biomix_analysis.R
# Main analysis pipeline using modular BiomiX functions

# =======================
# User Parameters Manual load (Debugging)
# =======================
#
# library(vroom)
# args = as.list(c("Neutrophils","PAPS"))
# args[1] <-"RA"
# args[2] <-"CTRL"
# args[3] <-"C:/Users/crist/Desktop/BiomiX2.5"
# 
# directory <- args[[3]]
# iterations = 1
# selection_samples = "NO"
# Cell_type = "Plasma"
# i = 2
# ANNOTATION = "Annotated"
# DIR_METADATA <- readLines("C:/Users/crist/Desktop/BiomiX2.5/directory.txt")
# STATISTICS = "YES"
# renv::load(paste(directory,"_INSTALL",sep="/"))

# =======================
# Load Functions & Setup
# =======================

library(vroom)
library(dplyr)
library(cmmr)
library(stringr)
library(rlist)
library(tibble)
library(readxl)


setwd(directory)



#LOADING ARGUMENTS AND VARIABLES REQUIRED IN THE CODE


library(jsonlite)
library(tidyverse)
setwd(directory)
print(directory)
combined_json <- jsonlite::fromJSON(txt = "COMBINED_COMMANDS.json")


COMMAND <- combined_json[["COMMANDS"]]
COMMAND_MOFA <- combined_json[["COMMANDS_MOFA"]]
COMMAND_ADVANCED <- combined_json[["COMMANDS_ADVANCED"]]
DIR_METADATA <- combined_json[["DIRECTORY_INFO"]][["METADATA_DIR"]]
DIR_METADATA_output <- combined_json[["DIRECTORY_INFO"]][["OUTPUT_DIR"]]
Heatmap_genes <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_CLUSTERING_OPTIONS[3])
LogFC <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS[1])
padju <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS[2])
ANNOTATION <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL[1]
ANNOTATION_TYPE <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL[2]
MS2_databases = c("HMDB","MONA","MASSBANK")  


directory2 <- paste(directory,"/Metabolomics/",sep="")
setwd(directory2)
source("BiomiX_DMA_functions.r")








#ORDER THE MS2_DATABASES PRIORITY BASED ON THE USER INDICATION IN THE ADVANCED OPTIONS.
Input_ms2<- COMMAND_ADVANCED[3,8]
Input_ms2 <- unlist(strsplit(as.character(Input_ms2), "/"))
Input_ms2<-as.numeric(substr(Input_ms2, 1, 1))
if(sum(is.na(Input_ms2)) > 0){
        MS2_databases<-MS2_databases[!is.na(Input_ms2)]
        Input_ms2<-Input_ms2[!is.na(Input_ms2)]
}
MS2_databases <-MS2_databases[order(Input_ms2)]


directory2 <- paste(directory,"/Metabolomics/INPUT",sep="")
setwd(directory2)

print(i)
print(COMMAND$LABEL[i])




### LOADING ANNOTATION AND PEAKS MATRIX
if (grepl("\\.xlsx$|\\.xls$", COMMAND$DIRECTORIES[i])) {
        matrix <- read_excel(COMMAND$DIRECTORIES[i])
        print("Matrix Excel File read successfully!")
}else{
        matrix <- vroom(COMMAND$DIRECTORIES[i],delim="\t",col_names = TRUE, comment = "#")
}


#Load Mart metabolome data
Mart_metabolome<-vroom(paste(directory,"/Integration/x_BiomiX_DATABASE/Metabolome_reference_ID.tsv", sep=""), delim = "\t")
if (ANNOTATION == "Annotated" | ANNOTATION_TYPE == "compound_name") {
  matrix <- matrix %>%
    left_join(Mart_metabolome %>% select(name, HMDB),
              by = c("ID" = "name")) %>%
    mutate(ID = if_else(!is.na(HMDB), HMDB, ID)) %>%
    select(-HMDB)
  matrix <- matrix[grep("^HMDB.*",matrix$ID),]
}

sam <- colnames(matrix)
pea <-matrix$ID
matrix <- t(matrix[,-1])
matrix <- add_column(as.data.frame(matrix), sam[-1], .after = 0) 
colnames(matrix) <-c("ID",pea)



#LOAD METADATA

if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
  Metadata_total <- read_excel(DIR_METADATA)
  print("Metadata Excel File read successfully!")
}else{
  Metadata_total <- vroom(DIR_METADATA, delim = "\t", col_names = TRUE)}






### UPLOAD ANNOTATION IF THE USER SELECTED THE MS1 OR MS2 ANNOTATION
if(ANNOTATION != "Annotated"){
ANNOTATION = COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL[1]

load_annotation(ANNOTATION, i, COMMAND_ADVANCED)



#create directory
dir.create(path = paste(directory,"/Integration/INPUT/", "Metabolomics_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
directory2 <- paste(directory, "/Integration/INPUT/", "Metabolomics_", Cell_type,  "_",args[1],"_vs_", args[2], sep ="")

setwd(directory2)


### HMDB METABOLITE FILTERING  
#UPLOAD PLASMA, URINE OR OTHERS ANNOTATIONS FROM HMDB AND FILTERING 
#IT WORKS IF THE LABEL IS NAMES AS THE TISSUES/SAMPLE WHERE THE DATABASE IS AVAILABLE 


# Example calls (equivalent to your original code)
load_biofluid_metabolites(COMMAND$LABEL[i], "Plasma", "serum_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites.csv?action=index&blood=1&c=hmdb_id&controller=metabolites&d=up&detected=1&expected=1&filter=true&predicted=1&quantified=1&utf8=%E2%9C%93",
                          "serum_metabolite", directory2)

load_biofluid_metabolites(COMMAND$LABEL[i], "Urine", "urine_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites.csv?action=index&c=hmdb_id&controller=metabolites&d=up&detected=1&expected=1&filter=true&predicted=1&quantified=1&urine=1&utf8=%E2%9C%93",
                          "urine_metabolite", directory2)

load_biofluid_metabolites(COMMAND$LABEL[i], "Saliva", "saliva_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&saliva=1&filter=true",
                          "saliva_metabolite", directory2)

load_biofluid_metabolites(COMMAND$LABEL[i], "Cerebrospinal Fluid", "cerebrospinal_fluid_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&csf=1&filter=true",
                          "CSF_metabolite", directory2)

load_biofluid_metabolites(COMMAND$LABEL[i], "Feces", "feces_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&feces=1&filter=true",
                          "feces_metabolite", directory2)

load_biofluid_metabolites(COMMAND$LABEL[i], "Sweat", "sweat_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&sweat=1&filter=true",
                          "sweat_metabolite", directory2)

load_biofluid_metabolites(COMMAND$LABEL[i], "Breast milk", "breast_milk_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&breast_milk=1&filter=true",
                          "breast_milk_metabolite", directory2)

load_biofluid_metabolites(COMMAND$LABEL[i], "Bile", "bile_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&bile=1&filter=true",
                          "bile_metabolite", directory2)

load_biofluid_metabolites(COMMAND$LABEL[i], "Amniotic Fluid", "AF_metabolite_annotated.tsv",
                          "https://hmdb.ca/metabolites?utf8=%E2%9C%93&filter=true&quantified=1&detected=1&expected=1&predicted=1&amniotic_fluid=1&filter=true",
                          "AF_metabolite", directory2)


}

#SAMPLES SELECTUON BASED ON THEIR LABEL NAME (BY REGEX) -
# IT CORRESPONDS TO -> SELECTION OPTION ON INTERFACE
process_sample_selection(selection_samples, Cell_type, Metadata_total, matrix)


Metadata <- Metadata_Bcell
Metadata_individual=NULL
Metadata_reads=NULL
Metadata_Bcell=NULL
Metadata_tot = NULL
Metadata_total = NULL

print("Same samples in matrix and metadata?")
print(all(Metadata$ID == matrix[,1]))
outs <- Metadata$CONDITION  == args[2]|Metadata$CONDITION  == args[1] | Metadata$CONDITION  == "QC"
Metadata<- Metadata[outs,]
matrix<- matrix[outs,]


#FILTERING SAMPLES BASED ON METADATA CRITERIA SELECTED. 
#DEFINED IN ADVANCED OPTION, METADATA SECTION

filter_samples_by_metadata(COMMAND_ADVANCED, Metadata, matrix)



###
matrix<- as.data.frame(matrix)
rownames(matrix) <- matrix$ID
###


#This part of code is to refill the column names having an NA name and rename the duplicated ones
na_cols <- which(is.na(colnames(matrix)))
new_names <- paste("UNKNOWN",1:sum(is.na(colnames(matrix))),sep="") 
colnames(matrix)[na_cols] <- new_names
colnames(matrix) <- make.unique(names(matrix), sep = "_")


#QC VISUALISATION FUNCTION
if(COMMAND$PREVIEW[i] == "YES"){
browser_analysis <- readLines(paste(directory, '/_INSTALL/CHOISE_BROWSER_pre-view', sep = ""), n = 1)
run_qc_preview(COMMAND, i, matrix, Metadata, directory)
} else {
  print("no QC pre-visualization")
}


#SAVING THE INPUT FOR THE MOFA ANALYSIS

matrixs <- matrix 
matrixs <- add_column(matrixs, Metadata$CONDITION, .after = 1)
colnames(matrixs)[2] <- "CONDITION"
matrixs[,3:ncol(matrixs)]<-apply(matrixs[,3:ncol(matrixs)],2,as.numeric)
#matrixs[matrixs == 0] <- 1


dir.create(path = paste(directory,"/Integration/INPUT/", "Metabolomics_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")

directory3 <- paste(directory, "/Integration/INPUT/", "Metabolomics_", Cell_type,  "_",args[1],"_vs_", args[2], sep ="")

write.table(matrixs,paste(directory3,"/Metabolomics_",Cell_type, "_MOFA.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")


if (STATISTICS == "YES"){


# #### FUNCTION DEFINITION FOR STATISTICAL TEST BETWEEN CONTROL AND TESTED SAMPLES----



# #### DMS DIFFERENTIAL METABOLITE SIGNALS ANALYSIS----

CON <-DMS_REFERENCE(args[2])        
TEST <- DMS_TESTED(args[1])        
TEST$NAME <- row.names(TEST)
gc()

#CHECK IF THE PEAKS ARE ALREADY ANNOTATED OR NOT, IF NOT IT WILL COMPARE THE PEAK
#SIGNALS DIRECTLY, OTHERWISE IT WILL START THE ANNOTATION AT FIRST.
if(ANNOTATION == "Annotated"){

###  DGE ANALYSIS ###  
run_dge_analysis_annotated(TEST, matrix, Metadata, Cell_type, args, directory, padju, LogFC)

###  HEATMAP ###
generate_heatmap_annotated(ALTI, BASSI, matrix, Metadata, Heatmap_genes, Cell_type, args, directory2)

query_id <- total[which(total$p_val < 0.05),] #REPLACE IT WITH P.ADJ

##### MetPath  #####
run_metpath_pipeline_annotated(total, COMMAND_ADVANCED, Cell_type, args, directory2, hmdb_pathway, kegg_hsa_pathway)

##### MetaboAnalistR  #####
run_metaboanalyst_pipeline_annotated(query_id, COMMAND_ADVANCED, Cell_type, args, directory, directory2)        
        
} else {
        
        if(ANNOTATION == "MS2"){  
          
          
          #Generation parameters for tidymass annotation multidataset
          run_msms_filtering(COMMAND_ADVANCED, directory, MS2_databases)
                
                param2 <- param
                
                #RENAME DATABASES TO UPPERCASE
                for (it in 1:length(param)){
                        if( MS2_databases[it] == toupper(param[[1]]$database@database.info$Source)){
                                param2[[it]] <- param[[1]]
                        }else if(MS2_databases[it] == toupper(param[[2]]$database@database.info$Source)){
                                param2[[it]] <- param[[2]]
                        }else{
                                param2[[it]] <- param[[3]]  
                        }}
                
                
                param <- param2
                
                gc()
                setwd(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY[1])
                print(annotation)
                annotation2 <- annotation[,-4]
                colnames(annotation2) <- c("\"name\"","\"mz\"","\"rt\"")
                write.table(annotation2,"TEMP.csv" ,quote = FALSE, row.names = F, sep = ",")
                
                #FIND THE LIST OF mzML file belonging to each matrix and add all of them in a vector
                files <- grep("MS2*",list.files(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY[1]),value=TRUE)
                files_MS2<-grep(paste("MS2_",COMMAND$LABEL[i],"*",sep=""),files,value=TRUE,ignore.case =TRUE)
                
                Sys.time()
                setwd(directory3)
                

                
                #MULTI DATASET ANNOTATION BY TIDYMASS
                print(files_MS2)
                annotate_result5 <- 
                        identify_metabolite_all(ms1.data = "TEMP.csv", 
                                                ms2.data = files_MS2, #USING ONLY ONE FILE
                                                parameter.list = param,
                                                path = path)
                
                Sys.time()
                
                setwd(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY[1])
                
                if (file.exists("TEMP.csv")) {
                        file.remove("TEMP.csv")
                }
                
                ## RESTORE ####
                if (dir.exists("intermediate_data")) {
                        unlink("intermediate_data", recursive = TRUE)
                }
                
                ## RESTORE ####
                if (dir.exists("Result")) {
                        unlink("Result", recursive = TRUE)
                }
                
                
                
                
                
                #This part of the script is made to select all the peaks found in all the
                #databases. The final file (annot) will be used to replace the metabolites
                #with a level of annotation of 3 or 4.
                
                build_msms_annotation(directory, Cell_type, args, param, annotate_result5)
                
                #Genetartion MS2 spectra for each peak and the metabolites identified in the reference databases
                generate_ms2_spectra_plots(directory, directory2, Cell_type, args, MS2_databases, annotate_result5)
                  
                
        }else{
                
                #This is the starting of the MS1 pipeline and it is shared with the MS2
                dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                dir.create(path = paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                directory2 <- paste(directory,"/Metabolomics/OUTPUT/", Cell_type, "_", args[1], "_vs_", args[2], sep ="")
                
                
                setwd(directory2)        
        }
        
        gc()
        
        total <- merge(TEST,annotation,by.x = "NAME", by.y = "name",all.x = TRUE)
        
        total$padj = p.adjust(total$p_val, method = "fdr") 
        
        total_min <- total[,c("NAME", "log2FC", "p_val", "padj")]
        total_min <- total_min %>% arrange(padj)
        write.table(total_min,paste(directory2,"/",Cell_type,"_",args[1],"_vs_",args[2],"_peak_statistics.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")
        
        
        
        library(stringr)
        
        
        #Retrieve the annotation using the advanced research in ceu database
        
        if(ANNOTATION == "MS1"){ 
                adduct_list <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2[1]
                mode_ions <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1[1]
                adduct_list_2 <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2[2]
                tolerance_list <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1[3]
                list_database_ms1 <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2[3]}
        
        if(ANNOTATION == "MS2"){ 
                adduct_list <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2[1]
                mode_ions <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3[1]
                adduct_list_2 <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2[2]
                tolerance_list <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4[2]
                list_database_ms1 <- COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4[1]}
        
        #==============================================================================
        if(mode_ions == "positive" | mode_ions == "negative"){
                
                Positive_adduct <- as.array(adduct_list)
                Positive_adduct <- strsplit(Positive_adduct, "/")
                Positive_adduct <- Positive_adduct[[1]]
                
                Negative_adduct <- as.array(adduct_list_2)
                Negative_adduct <- strsplit(Negative_adduct, "/")
                Negative_adduct <- Negative_adduct[[1]]
                
                if (is.na(Negative_adduct)){
                        ducts <- as.list(Positive_adduct)
                        #ducts <- Positive_adduct %>% str_c(collapse = ",")
                        #ducts <- paste("[",ducts,"]",sep = "")
                } else if (is.na(Positive_adduct)){
                        ducts <- as.list(Negative_adduct)
                        #ducts <- Positive_adduct  %>% str_c(collapse = ",") 
                        #ducts <- paste("[",ducts,"]",sep = "")
                } else{
                        #ducts <- c(Positive_adduct,Negative_adduct)
                        print("ERROR YOU SELECTED POSITIVE AND NEGATIVE ADDUCTS AT THE SAME TIME!!")
                }
                
                masses_mode <- as.array("mz")
                Ion_mode <- as.array(mode_ions)
                
        } else if(mode_ions == "neutral"){
                Ion_mode <- as.array("neutral")
                masses_mode <- as.array("neutral")
                ducts <- "[\"all\"]"
        }
        
        tolerance <- as.array(tolerance_list)
        
        
        databases <- as.array(list_database_ms1)
        databases <- strsplit(databases, "/")
        databases <- databases[[1]]
        databases <- as.list(databases)
        #databases <- databases %>% str_c('"', ., '"') %>% str_c(collapse = ",") 
        #databases<-paste("[",databases,"]",sep = "")
        
        
        
        #MASS <- as.array(round(as.numeric(total$`m/z`),4))
        MASS <- round(as.numeric(total$`m/z`),4)
        #RT <- as.array(total$RT_min)
        RT <- as.numeric(total$RT_min)
        
        # MASS <- MASS[0:300]
        # RT <- RT[0:300]
        # RT <- rep(1,300)
        
        
        cat(MASS, fill = getOption("width"), sep = ",")
        # xx<-paste(MASS, collapse = " ")
        # xx<-gsub(" ", ",", xx, fixed=TRUE)
        # MASS<-paste("[",xx,"]",sep = "")
        
        cat(RT, fill = getOption("width"), sep = ",")
        # xx<-paste(RT, collapse = " ")
        # xx<-gsub(" ", ",", xx, fixed=TRUE)
        # RT<-paste("[",xx,"]",sep = "")
        
        
        #THIS IF ELSE ALLOWS TO UPLOAD A PREVIOUS ANNOTATION MADE ON THE SAME DATASET.
        #IF THIS ANNOTATION IS NOT AVAILABLE OR IF IT IS THE FIRST ANALYSIS WILL AS TO
        #CEU MASS MEDIATOR TO DO IT.
        
        retrieve_advanced_batch_annotation(directory2, Cell_type, args, databases, masses_mode, Ion_mode, ducts, tolerance, MASS, RT)
        
        #===============================================================================
        ### RESTART FROM HERE
        
        
        
        
        #save the annotation to avoid to run again the annotation search
        
        #===================================UPLOAD SAVED FILE==============================================
        
        head(advanced_batch_df)
        str(advanced_batch_df)
        gc()
        
        as.numeric(total$`m/z`) %in% advanced_batch_df$experimental_mass
        
        #Match annotation file and peaks file (total) using the m/z to match
        
        advanced_batch_df$experimental_mass <- as.character(round(as.numeric(advanced_batch_df$experimental_mass),4))
        total$`m/z` <- as.character(round(as.numeric(total$`m/z`),4))
        total$`m/z`<-as.character(total$`m/z`)
        
        
        
        tat<-merge(total,advanced_batch_df, by.x="m/z", by.y="experimental_mass", all.x=TRUE)
        total$`m/z` %in% advanced_batch_df$experimental_mass
        
        
        x <- colnames(tat) %in% matrix$ID
        tat <- tat[,!x]
        #merge the metabolites found with the table containing the results
        
        
        print(i)
        print(COMMAND$LABEL[i])
        filtering= 0
        
        

        
        
        #### FILTER OF TISSUES SPECIFIC METABOLITES IN HMDB #####
        #IF THE ANNOTATION FOR SERUM OR URINE IS AVAILABLE IT IS USED
        #TO FILTER THE RESULTS OBTAINED USING CEU MASS MEDIATOR
        
        result <- NULL
        
        # Try each biospecimen
        result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, serum_metabolite, "Serum")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, serum_metabolite, "Plasma")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, urine_metabolite, "Urine")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, saliva_metabolite, "Saliva")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, CSF_metabolite, "Cerebrospinal Fluid")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, Feces_metabolite, "Feces")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, sweat_metabolite, "Sweat")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, breast_milk_metabolite, "Brest Milk")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, bile_metabolite, "Bile")
        if (is.null(result)) result <- filter_by_biospecimen(COMMAND$LABEL[i], tat, AF_metabolite, "Amniotic Fluid")
        
        # Fallback: no filtering applied
        if (is.null(result)) {
          print("NO BIOSPECIMENS filtering")
          total <- tat
          filtering <- "NO"
          colnames(total)[which(colnames(total) == "NAME")] <- "NAME.x"
        } else {
          total <- result$total
          filtering <- result$filtering
        }
        
        
        
        totalshOK <- total %>% filter(padj <padju)
        total <- total %>% arrange(padj)
        print(colnames(total))
        
        # x <- colnames(total) %in% matrix$ID
        # total <- total[,!x]
        
        total$MS2_annot <- "NA"
        total$annot_level <- "Level_3"
        
        to_eliminate=NULL
        to_eliminate2=NULL
        
        if(ANNOTATION == "MS2"){  
                
                
          result <- replace_ms1_with_ms2(annot, total, to_eliminate2)
          total <- result$total
          to_eliminate2 <- result$to_eliminate2
          
                
        }
        
        
        write.table(total,paste(directory2,"/",Cell_type,"_",args[1],"_vs_",args[2],"_results.tsv", sep = ""),quote = FALSE, row.names = F, sep = "\t")
        
        
        
        
        pdf(file=paste("plot_DMA_", args[1],"_",args[2],"_", Cell_type,".pdf",sep=""))
        
        dge_results <- perform_dge_analysis_MS1_MS2(total, padju, LogFC, Cell_type, args, Metadata, directory2)
        ALTI <- dge_results$ALTI
        BASSI <- dge_results$BASSI
        
        generate_heatmap_MS1_MS2(ALTI, BASSI, matrix, Heatmap_genes, Metadata, Cell_type, args, directory2)
        
          
          
        ##### MetPath  #####
        library(metpath)
        library(tidyverse)
        library(dplyr)
        
        dir.create(path = paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        setwd(paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep =""))
        
        saved<-total
        #Filtering of peak with one or max 2 annotation
        peak_repeated <-total %>% count(NAME.x)
        total<- total[!total$NAME.x %in% peak_repeated$NAME.x[peak_repeated$n > 2],]
        total <- total %>% arrange(NAME.x)
        
        #Eliminate level 3 peaks with level 2 annotation found
        peak_repeated<-total$NAME.x[total$annot_level == "Level_2"]
        outs<-total$NAME.x %in% peak_repeated & total$annot_level == "Level_3"
        total <-total[!outs,]
        
        
        
        #Filtering annotation couples where both or a single one contains the annotation
        #Selection of best annotation in peak with 2 annotation based on minimun ppm error
        iter<-which(duplicated(total$NAME.x) & !total$annot_level == "Level_2")
        to_elim = NULL
        iterat=NULL
        for (u in iter){
                n1<-u
                s1<-total$error_ppm[u]
                n2<-u-1
                s2<-total$error_ppm[u-1]
                
                #Statement to treat missing ppm_errors values
                if (is.na(s1) | is.na(s2)){
                        if(is.na(s1) & is.na(s2)){
                                to_elim<-append(to_elim, n2)                                
                        }else if(is.na(s1)){
                                to_elim<-append(to_elim, n1)
                        }else{#(is.na(s2))
                                to_elim<-append(to_elim, n2)
                        }
                        
                }else{
                
                
                if(as.numeric(s1) == as.numeric(s2)){
                        #The results with the same ppm error are kept for the metpath analysis
                        }else{
                        
                elim <- max(as.numeric(s1), as.numeric(s2))
                if(elim == as.numeric(s1)){
                        print = as.numeric(s1)
                        to_elim<-append(to_elim, n1)
                }
                if(elim == as.numeric(s2)){
                        to_elim<-append(to_elim, n2)
                }
                
                        }
                        }
                
        }
        
        if(length(to_elim) != 0){
                total <- total[-to_elim, ]}
        
        
        
        
        run_metpath_pipeline_MS1_MS2(total,
                                 saved,
                                 directory2,
                                 Cell_type,
                                 args,
                                 hmdb_pathway,
                                 kegg_hsa_pathway,
                                 COMMAND_ADVANCED)
        
        
        
        query_id <- total[which(total$p_val < 0.05),] #REPLACE IT WITH P.ADJ
        query_id_saved <- saved[which(saved$p_val < 0.05),]
        
        #SAVE PATHWAY ANALYSIS RESULTS
        if (nrow(query_id) != 0){
                
                dir.create(path = paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                setwd(paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], sep =""))
                
                
                pathway_class_HMDB = 
                        metpath::pathway_class(hmdb_pathway)
                pathway_class_KEGG = 
                        metpath::pathway_class(kegg_hsa_pathway)
                
                
                pdf(file=paste("Pathway_analysis_HMDB_", args[1],"_",args[2],"_", Cell_type,".pdf", sep=""))
                
                
                gc()
                remain_idx = which(unlist(pathway_class_HMDB) == "Metabolic;primary_pathway")
                hmdb_pathway = hmdb_pathway[remain_idx]
                hmdb_pathway
                
                
                result = 
                        enrich_hmdb(query_id = unique(query_id$HMDB), 
                                    query_type = "compound", 
                                    id_type = "HMDB",
                                    pathway_database = hmdb_pathway,
                                    only_primary_pathway = TRUE,
                                    p_cutoff = 0.05, 
                                    p_adjust_method = "BH", 
                                    threads = as.numeric(COMMAND_ADVANCED[3,3]))
                
                result
                
                
                if (length(result) != 0){
                        
                        x<-enrich_bar_plot(
                                object = result,
                                x_axis = "p_value_adjust",
                                cutoff = 1.1,
                                top = 10
                        )
                        
                        print(x)
                        
                        x <-enrich_scatter_plot(object = result)
                        
                        print(x)
                        
                        write.table(result@result,paste(directory2,"/Pathway_analysis/", "HMDB_", Cell_type, "_",args[1],"_vs_", args[2], "/HMDB_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
                        
                        dev.off()
                        
                }
                
                gc()
                dir.create(path = paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
                setwd(paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], sep =""))
                
                
                pdf(file=paste("Pathway_analysis_KEGG_", args[1],"_",args[2],"_", Cell_type,".pdf", sep=""))
                
                
                head(pathway_class_KEGG)
                remain_idx =
                        pathway_class_KEGG %>%
                        unlist() %>%
                        stringr::str_detect("Disease") %>%
                        `!`() %>%
                        which()
                
                remain_idx
                
                pathway_database =
                        kegg_hsa_pathway[remain_idx]
                
                pathway_database
                
                
                result = 
                        enrich_kegg(query_id = unique(query_id$KEGG), 
                                    query_type = "compound", 
                                    id_type = "KEGG",
                                    pathway_database = pathway_database, 
                                    p_cutoff = 0.05, 
                                    p_adjust_method = "BH", 
                                    threads = as.numeric(COMMAND_ADVANCED[3,3]))
                
                
                result
                
                if (length(result) != 0){
                        
                        x <-enrich_bar_plot(
                                object = result,
                                x_axis = "p_value_adjust",
                                cutoff = 1.1,
                                top = 10
                        )
                        print(x)
                        
                        x<-enrich_scatter_plot(object = result)
                        print(x)
                        
                        write.table(result@result,paste(directory2,"/Pathway_analysis/", "KEGG_", Cell_type, "_",args[1],"_vs_", args[2], "/KEGG_table_results.tsv", sep =""),quote = FALSE, row.names = F, sep = "\t")
                        
                        dev.off()
                        
                }
                
        }else{
                print("NO STATISTICAL SIGNIFICANT PEAK FOR METABOLITE PATHWAY ANALYSIS")
        }
        
        
        ##### MetaboAnalistR  #####
        
        
        run_metaboanalyst_pipeline_MS1_MS2(query_id,
                                            query_id_saved,
                                            total,
                                            Cell_type,
                                            args,
                                            directory,
                                            directory2
                                            )
        
                
}

}else{
        print("No statistical analysis and annotation")
}


advanced_batch_df = NULL
param=NULL
param1=NULL
annotation = NULL
annotation2 =NULL
massbank_database0.0.3 =NULL
mona_database0.0.3= NULL
pathway_class_KEGG =NULL
pathway_class_HMDB =NULL
hmdb_database0.0.3=NULL
tat= NULL
NO=NULL
serum_metabolite =NULL
annotate_result5=NULL
results=NULL
gc()

