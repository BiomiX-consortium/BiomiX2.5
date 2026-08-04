# Code to integrate data using SNF and evaluate the clustering results
# Author: Jessica Gliozzo

# Load libraries----
library("dplyr");
library("pheatmap");
library("igraph");
library("visNetwork");
library("htmlwidgets");
library("fpc"); # internal validation indices
library("ggplot2");
library("tidyr");    
#library("mclustcomp"); # external validation indices
library("aricode"); # ARI, NMI, AMI, NVI
# library("sabre"); # V-measure
library("ggalluvial");
library("survival"); # log-rank test and KM
library("survminer"); # plot KM
library("purrr"); # logging 
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
library("SNFtool");
library("RColorBrewer");
library("vroom");
library("rlist");
library("readxl");
library("jsonlite")
library("tidyverse")

#### STANDALONE ENTRY POINT (only when NOT source()'d) ----
#
# Invocation (Docker mode, via runner_SNF in BiomiX_BETA.r):
#   Rscript SNF_int.R <DIRECTORY> <CELL_TYPE> <GROUP_1> <GROUP_2> <SHARED_DIR> <ITERATIONS>
#
# Local mode: source()'d from BiomiX_BETA.r, which already has its own
# `directory` and 3-element `args` (group1, group2, directory) in scope.

if (!exists("Cell_type")) {
  
  cli_args <- commandArgs(trailingOnly = TRUE)
  if (length(cli_args) < 6) {
    stop("Usage: Rscript SNF_int.R <DIRECTORY> <CELL_TYPE> <GROUP_1> <GROUP_2> <SHARED_DIR> <ITERATIONS>")
  }
  
  directory  <- cli_args[1]
  args       <- c(cli_args[3], cli_args[4])
  shared_dir <- cli_args[5]
  Cell_type  <- cli_args[2]
  
  #renv::load(paste(directory, "_INSTALL", sep = "/"))
  
  # site_lib <- "/usr/local/lib/R/site-library"
  # if (dir.exists(site_lib) && !(site_lib %in% .libPaths())) {
  #   .libPaths(c(.libPaths(), site_lib))
  # }
  
  combined_json <- jsonlite::fromJSON(txt = file.path(shared_dir, "COMBINED_COMMANDS.json"))
  
} else {
  # Local/sourced case: `directory` and 3-element `args` already exist from
  # BiomiX_BETA.r. Only shared_dir needs a fallback if not already present.
  if (!exists("shared_dir")) shared_dir <- directory
  directory <- unlist(args[[3]])
  setwd(directory)
  combined_json <- jsonlite::fromJSON(txt = file.path(shared_dir, "COMBINED_COMMANDS.json"))
}

int_method <- "SNF"

MART <- vroom(paste(directory,"/Integration/x_BiomiX_DATABASE/mart_export_37.txt",sep=""), delim = ",")
myList <- list()

COMMAND <- combined_json[["COMMANDS"]]
COMMAND_MOFA <- combined_json[["COMMANDS_MOFA"]]
COMMAND_ADVANCED <- combined_json[["COMMANDS_ADVANCED"]]

Max_features_SNF <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[3])

K.snf <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[1])
sigma <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[2])
t <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[3])
nc <- c(2:as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[1]))

apply_scaling_SNF = as.logical(COMMAND_ADVANCED$ADVANCED_OPTION_NEMO_OPTIONS[2])

if (nc[length(nc)] < 2) {
  stop("Maximum number of clusters must be at least 2.")
}

top_feat <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[2])
enrich_vars <- trimws(strsplit(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[1]), "/")[[1]])
surv_vars <- trimws(strsplit(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[2]), "/")[[1]])
DIR_METADATA <- combined_json[["DIRECTORY_INFO"]][["METADATA_DIR"]]

gt.clust_name <- trimws(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[3]))

directory2 <- paste(directory,"/Integration",sep="")

if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
  METADATA <- read_excel(DIR_METADATA)
  print("Metadata Excel File read successfully!")
}else{
  METADATA <-vroom(DIR_METADATA , delim = "\t", col_names = TRUE)}

set.seed(123)
setwd(directory2)

# Load functions----
source("PSN_utils.R")


int_method <- "SNF"

MART <- vroom(paste(directory,"/Integration/x_BiomiX_DATABASE/mart_export_37.txt",sep=""), delim = ",")
myList <- list()

#combined_json file loaded from the BiomiX_BETA.r script (master script) 
#and loading the COMBINED_COMMANDS.json file (master commands).

COMMAND <- combined_json[["COMMANDS"]]
COMMAND_MOFA <- combined_json[["COMMANDS_MOFA"]]
COMMAND_ADVANCED <- combined_json[["COMMANDS_ADVANCED"]]

#Max_features <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[3])
Max_features_SNF <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[3])

K.snf <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[1]) # number of neighbors in KNN
sigma <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[2]) # variance for affinityMatrix
t <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[3]) # number of iterations for SNF
nc <- c(2:as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[1])) # Max number of cluster to test

#To be included in the interface:
apply_scaling_SNF = as.logical(COMMAND_ADVANCED$ADVANCED_OPTION_NEMO_OPTIONS[2]) # JG: it could be set as "TRUE"/"FALSE"

if (nc[length(nc)] < 2) {
    stop("Maximum number of clusters must be at least 2.")
}

top_feat <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[2]) #Number variable in the heatmap to visualize
# Variable of interest for enrichment and survival analysis
enrich_vars <- trimws(strsplit(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[1]), "/")[[1]])
surv_vars <- trimws(strsplit(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[2]), "/")[[1]])
DIR_METADATA <- combined_json[["DIRECTORY_INFO"]][["METADATA_DIR"]]

# Ground truth clustering name (if available, otherwise NULL)
gt.clust_name <- trimws(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[3]))

# MANUAL INPUT JESS
# DIR_METADATA <- "~/Documents/biomix_project/BiomiX2.5/Metadata/EGAS00001001746_metadata_CLL.tsv"
# COMMAND$DATA_TYPE <- c("Transcriptomics", "Methylomics", "Undefined", "Undefined", "Undefined", "Undefined", "X")
# COMMAND$INTEGRATION <- c("YES","YES", "NO", "NO","NO","NO") #CLL (transcriptomics + methylomics)
# COMMAND$LABEL <- c("RNA",  "Methylomics", NA, NA, NA, NA)

directory2 <- paste(shared_dir,"/Integration",sep="")

if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
  METADATA <- read_excel(DIR_METADATA)
  print("Metadata Excel File read successfully!")
}else{
  METADATA <-vroom(DIR_METADATA , delim = "\t", col_names = TRUE)}


# Clean workspace, set seed, set working directory----
# rm(list=ls())
set.seed(123)
setwd(directory2)




# Load functions----
#source("PSN_utils.R");
# # Collect hyperparameters for the analysis (variables and values used defined)----
# int_method = "SNF" # integration method (can be "SNF" or "NEMO")
# 
# K.snf = 10 # number of neighbors in KNN
# sigma = 0.5 # variance for affinityMatrix
# t = 15 # number of iterations for SNF
# 
# # Variable of interest for enrichment and survival analysis
# enrich_vars <- c("GENDER", "AGE")
# surv_vars <- c("OS.time", "OS.event", "DFS.time", "DFS.event")
# 
# # Ground truth clustering name (if available, otherwise NULL)
# gt.clust_name <- "CONDITION"
# 
# # Range of clusters to test
# nc=2:8
# 
# # Top features to plot in heatmap
# top_feat = 20

# MANUAL INPUT JESS
# directory="~/Documents/biomix_project/BiomiX2.5"
# DIR_METADATA <- "~/Documents/biomix_project/BiomiX2.5/Metadata/EGAS00001001746_metadata_CLL.tsv"
# COMMAND$DATA_TYPE <- c("Transcriptomics", "Methylomics", "Undefined", "Undefined", "Undefined", "Undefined", "X")
# COMMAND$INTEGRATION <- c("YES","YES", "NO", "NO","NO","NO", "NO") #CLL (transcriptomics + methylomics)
# COMMAND$LABEL <- c("RNA",  "METHY", NA, NA, NA, NA, NA)

##### REARRANGEMENT INPUT1 DATA ----

myList <- list()
names_X <- c()
for (i in 1:length(COMMAND$INTEGRATION)){
  
  
  #i <- 3 #N_input
  type <- COMMAND$DATA_TYPE[i]
  

if(COMMAND$INTEGRATION[i] == "YES"){
  if(COMMAND$DATA_TYPE[i] == "Metabolomics"){
    
    #directory2 <- paste(directory,"/Metabolomics",sep="")
    directory2 <- paste(shared_dir,"/Integration/INPUT/", "Metabolomics_", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
    serum_metabolomics <- vroom(paste(directory2,"/Metabolomics_",COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
    directory2 <- paste(shared_dir,"/Metabolomics/OUTPUT/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
    serum_annotation <- vroom( paste(directory2,"/",COMMAND$LABEL[i],"_",args[1],"_vs_",args[2],"_results.tsv", sep = ""), delim = "\t")
    INPUTX<-Metabolomics_processing(serum_annotation,serum_metabolomics)
    assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
    myList <- list.append(myList,INPUTX[[1]])
    names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])    
  }
  
  if(COMMAND$DATA_TYPE[i] == "Transcriptomics"){
    
    print(args[1])
    directory2 <- paste(shared_dir,"/Integration/INPUT/", COMMAND$LABEL[i],"_",args[1],"_vs_", args[2], sep ="")
    Wholeblood_RNAseq <-  vroom(paste(directory2, "/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], "_normalized_vst_variance.tsv",sep = ""), delim = "\t") #read normalization only
    Wholeblood_metadata <-  vroom(paste(directory2, "/","/Metadata_",COMMAND$LABEL[i], "_", args[1],".tsv",sep = ""), delim = "\t")
    INPUTX<-Transcriptomics_processing(Wholeblood_metadata,Wholeblood_RNAseq, Max_features_SNF)
    assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
    myList <- list.append(myList,INPUTX[[1]])
    names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])    
    
  }
  
  if(COMMAND$DATA_TYPE[i] == "Methylomics"){
    
    
    directory2 <- paste(shared_dir,"/Integration/INPUT/", "Methylome_",COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="") 
    Methylome_WB <-  vroom(paste(directory2, "/", COMMAND$LABEL[i], "_matrix_MOFA.tsv",sep = ""), delim = "\t") #read normalization only
    Methylome_metadata <-  vroom(paste(directory2, "/", COMMAND$LABEL[i],"_metadata_MOFA.tsv",sep = "") ,delim = "\t")
    directory2 <- paste(shared_dir,"/Methylomics/OUTPUT/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
    Methylome_annotation <- vroom(paste(directory2, "/", "DMP_", COMMAND$LABEL[i], "_Methylome_", args[1] ,"_vs_", args[2],".tsv",sep = ""), delim = "\t", col_names = TRUE)
    INPUTX<-Methylomics_processing(Methylome_annotation,Methylome_WB,Methylome_metadata, Max_features_SNF)
    assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
    myList <- list.append(myList,INPUTX[[1]])
    names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])    
  }
  
  if(COMMAND$DATA_TYPE[i] == "Undefined"){        
    directory2 <- paste(shared_dir,"/Integration/INPUT/", "Undefined_", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
    samples_undefined <- vroom(paste(directory2,"/Undefined_",COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
    INPUTX<-Undefined_processing(samples_undefined)
    assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
    myList <- list.append(myList,INPUTX[[1]])
    names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])
    
  }        
}
}

# Create directory to save all results ---
dir.create(path = paste(shared_dir,"/Integration/OUTPUT/",int_method, "_", args[1] ,"_vs_", args[2],sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")

directory2 <- paste(shared_dir,"/Integration/OUTPUT/",int_method, "_", args[1] ,"_vs_", args[2],sep="") 
setwd(directory2)

#Provide names to the list
names(myList) <- names_X


# Prepare data for SNF----
# Feature selection used only to test the code, not necessary in final script
data <- snf_nemo.preprocess(myList, METADATA, 
                            fsel=TRUE, Max_features = Max_features_SNF, 
                            apply_scaling=apply_scaling_SNF, 
                            int_method=int_method)

# Set ground truth clustering (if available) ----
# This should be selected by the user
if(!is.na(gt.clust_name) && gt.clust_name != "X"){
    # metadata are the same for all omics
    gt.clust <- as.factor(data$metadata[[gt.clust_name]]) # factor (ground truth clustering)
    names(gt.clust) <- data$metadata$ID
} else {
    gt.clust <- NULL # no ground truth clustering available
}

if (!gt.clust_name %in% colnames(data$metadata)) {
    stop("Ground-truth clustering variable not found in metadata: ", gt.clust_name)
}


# Compute similarity matrix for each modality----
sim_mat <- lapply(data$data, scaled.exp.euclidean, K=K.snf, sigma=sigma)

# save similarity matrices
dir.create("similarity_matrices", showWarnings = FALSE)
for(i in 1:length(sim_mat)){
    write.csv(sim_mat[[i]], file=file.path("similarity_matrices", 
              paste0(names(sim_mat)[i], "_similarity_matrix.csv")), 
              row.names=T)
}

# Check order of samples in metadata and similarity matrices
check_names(sim_mat, data$metadata)

# Integrate using SNF ----
W_int <- SNFtool::SNF(sim_mat, K=K.snf, t=t)


# Save integrated similarity matrix
write.csv(W_int, file=file.path("similarity_matrices", 
                                "integrated_similarity_matrix.csv"), row.names=T)


# Cluster optimization using eigengaps approach----

# Ground truth clustering (should be selected by user)
nc_estim <- estimate.nc(W_int, nc=nc, gt.clust=gt.clust, 
                        int_method=int_method)

# save optimal number of clusters
dir.create("clustering_metrics", showWarnings = FALSE)
write.csv(nc_estim, file="clustering_metrics/optimal_number_of_clusters.csv", 
          row.names=F)

# Spectral clustering----

# Apply clustering to all the number of clusters in nc_range
clustering <- lapply(nc_estim$nc_range, function(x) SNFtool::spectralClustering(W_int, K=x))
clustering <- as.data.frame(do.call(cbind, clustering))
rownames(clustering) = colnames(W_int)
colnames(clustering) = nc_estim$nc_range

# Save clustering results
dir.create("computed_clusterings", showWarnings = FALSE)
write.csv(clustering, file="computed_clusterings/clustering_results.csv", 
          row.names=T)

# Compute internal quality indices for clustering----
if(isSymmetric(W_int)){
  Wsym <- W_int  
} else {
    Wsym <- (W_int + t(W_int))/2 # enforce exact symmetry
}
smax <- max(Wsym[upper.tri(Wsym)])    # maximum value not considering diagonal
dist_mat <- 1 - (Wsym/smax)           # convert similarity to distance matrix
diag(dist_mat) <- 0                   # set diagonal to 0
int.val.idx <- apply(clustering, 2, function(x) fpc::cluster.stats(d=dist_mat, clustering=x, aggregateonly = T))
int.val.idx <- as.data.frame(do.call(cbind, int.val.idx))
to_remove <- which(rownames(int.val.idx) %in% c("g2", "g3", "corrected.rand", "vi"))
int.val.idx <- int.val.idx[-to_remove, , drop=F]


stability <- lapply(nc_estim$nc_range, function(x) fpc::clusterboot(dist_mat, B=100, distances=TRUE, 
                                                               bootmethod="subset", 
                                                               subtuning=round(0.8*nrow(dist_mat)), 
                                                               clustermethod=fpc_spectralClustering, 
                                                               seed=123, k=x, count=F))

# Extract subsetmean, subsetbrd and subsetrecover for each number of clusters in a dataframe
stability_results <- extract_subset_stats(stability, k_vec=nc_estim$nc_range)


# Visualize and save internal quality indices by number of clusters----
plot_int.val.idx(int.val.idx, nc_estim$nc_range)

write.csv(stability_results, file="clustering_metrics/stability_results.csv", row.names=F)

# Visualize unimodal and integrated similarity matrices----
# Possible types of visualizations are heatmap and graph

# create directory if does not exist
dir.create("heatmaps", showWarnings = FALSE)
dir.create("graphs", showWarnings = FALSE)

# make plots based on the optimal number of clusters obtained by eigengap approach
# or given as input (based on ground truth clustering).
# This is the "predicted/estimated" clustering.
pred.clust <- as.factor(clustering[, as.character(nc_estim$nc_estim)])
names(pred.clust) <- rownames(clustering)

# Plot heatmap
make_heatmap(c(sim_mat, list(W_int=W_int)), gt.clust=gt.clust, 
             pred.clust=pred.clust, norm="row", path="heatmaps")
# Plot graph
K_sparse <- rep(K.snf, length(sim_mat)+1)

make_graph(c(sim_mat, list(W_int=W_int)), gt.clust=gt.clust, 
           pred.clust=pred.clust, sparse_method="KNN", K=K_sparse, path="graphs")


# Compute external quality indices for clustering----
# External quality indices can be computed only if a ground truth clustering 
# is provided!!!
types <- c("ARI", "AMI", "NVI", "NMI")
if(!is.null(gt.clust)){
    y = as.numeric(gt.clust)
    ext.val.idx <- apply(clustering, 2,
                         function(x) {
                           cc <- aricode::clustComp(x, y);
                           res <- as.data.frame(setNames(
                             lapply(types, function(t) if (is.null(cc[[t]])) NA_real_ else cc[[t]]),
                             types
                           ));
                           rownames(res) <- NULL;
                           res})
    for(kn in names(ext.val.idx)){
        ext.val.idx[[kn]]$nc <- kn
    }
    ext.val.idx <- do.call(rbind, ext.val.idx)
    rownames(ext.val.idx) <- NULL
}


# Visualize and save external quality indices by number of clusters----
# External quality indices can be computed only if a ground truth clustering 
# is provided!!!
if(!is.null(gt.clust)){
    plot_ext.val.idx(ext.val.idx, nc_estim$nc_range, types=types)
}


# Make Alluvial plot to compare with ground truth clustering and variables----
# For now, I compare the optimal cluster to gt.clust (if present) and variables of interest (if present).

#print(enrich_vars)
valid_enrich_vars <- enrich_vars[enrich_vars != "X"]
if(!all(valid_enrich_vars %in% colnames(data$metadata))){
  stop("Some enrichment variables are not present in the metadata.")
}

if (length(valid_enrich_vars) > 0) {
  METADATA_alluvial <- as.data.frame(data$metadata[, valid_enrich_vars, drop=F])
  rownames(METADATA_alluvial) <- data$metadata$ID
  
  plot.alluvial(clustering[,as.character(nc_estim$nc_estim), drop=F], gt.clust, 
                METADATA_alluvial, 
                save_path="./clustering_metrics")

}

# Which are the most important features for integration?----
K.snf_vec <- rep(K.snf, length(data$data)) # K for each modality
feat_imp <- compute_feature_importance(data$data, pred.clust, K=K.snf_vec, sigma=sigma) 


dir.create("feature_importance", showWarnings = FALSE)
write.csv(feat_imp, file="feature_importance/feature_importance.csv", row.names=F)

# Heatmap of top_feat important features----
plot_imp_heatmap(data$data, feat_imp, pred.clust=pred.clust, gt.clust=gt.clust, 
                 top_feat=top_feat, save_path="./feature_importance")

# Compute variables enrichment and log-rank test----

# add fake survival data (just to test the function)
# set.seed(123)
# data$metadata$OS.time <- sample(1:100, nrow(data$metadata), replace=T)
# data$metadata$OS.event <- sample(0:1, nrow(data$metadata), replace=T)
# data$metadata$DFS.time <- sample(1:100, nrow(data$metadata), replace=T)
# data$metadata$DFS.event <- sample(0:1, nrow(data$metadata), replace=T)
# enrich_vars <- c("GENDER", "AGE")
# surv_vars <- c("OS.time", "OS.event", "DFS.time", "DFS.event")


dir.create("enrichment_survival_analysis", showWarnings = FALSE)


valid_surv_vars <- surv_vars[surv_vars != "X"]
# Keep only the variables with the following format .time / .event
valid_surv_vars <- valid_surv_vars[grepl("\\.time$|\\.event$", valid_surv_vars)]

if (length(valid_surv_vars) > 0 || length(valid_enrich_vars) > 0) {

  #quiet_enrich_surv_analysis <- quietly(enrich_surv_analysis) # Catch warnings
  enrich_surv_res <- enrich_surv_analysis(clustering, metadata=data$metadata, 
                                          valid_enrich_vars, valid_surv_vars, 
                                          file_path="enrichment_survival_analysis")
  
  enrich_surv_res_save <- rbind(enrich_surv_res$enrich_res, 
                                enrich_surv_res$surv_res)
  write.csv(enrich_surv_res_save, 
            file="enrichment_survival_analysis/enrichment_survival_results.csv", 
            row.names=T)
  
  write.csv(enrich_surv_res$warning_log, 
            file="enrichment_survival_analysis/enrichment_survival_results_warnings.csv", 
            row.names=F)
}

