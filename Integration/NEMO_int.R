# Code to integrate data using NEMO and evaluate the clustering results
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
library("NEMO");
library("RColorBrewer");
library("vroom");
library("rlist");
library("readxl");

# # #MANUAL INPUT
# args = as.list(c("BLymphocytes","SJS"))
# args[2] <-"CTRL"
# args[1] <-"RA"
# args[3] <- "C:/Users/crist/Desktop/BiomiX2.5"
# 
# directory <-args[3]
# 

# MANUAL INPUT JESS
# args = as.list(c("BLymphocytes","SLE"))
# args[1] <-"mutated"
# args[2] <-"unmutated"
# directory="~/Documents/biomix_project/BiomiX2.5"

int_method <- "NEMO"

MART <- vroom(paste(directory,"/Integration/x_BiomiX_DATABASE/mart_export_37.txt",sep=""), delim = ",")
myList <- list()

COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")
COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")

#Max_features <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[3])
Max_features_SNF <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[3])

K.nemo <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_NEMO_OPTIONS[1]) # number of neighbors in KNN
#sigma <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[2]) # variance for affinityMatrix
nc <- c(2:as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[1])) # Max number of cluster to test
if (nc[length(nc)] < 2) {
    stop("Maximum number of clusters must be at least 2.")
}

top_feat <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[2]) #Number variable in the heatmap to visualize
# Variable of interest for enrichment and survival analysis
enrich_vars <- trimws(strsplit(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[1]), "/")[[1]])
surv_vars <- trimws(strsplit(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[2]), "/")[[1]])
DIR_METADATA <- readLines(paste(directory,"directory.txt",sep="/"))

# Ground truth clustering name (if available, otherwise NULL)
gt.clust_name <- trimws(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[3]))


# MANUAL INPUT JESS
# DIR_METADATA <- "~/Documents/biomix_project/BiomiX2.5/Metadata/EGAS00001001746_metadata_CLL.tsv"
# COMMAND$DATA_TYPE <- c("Transcriptomics", "Methylomics", "Undefined", "Undefined", "Undefined", "Undefined", "X")
# COMMAND$INTEGRATION <- c("YES","YES", "NO", "NO","NO","NO", "NO") #CLL (transcriptomics + methylomics)
# COMMAND$LABEL <- c("RNA",  "METHY", NA, NA, NA, NA, NA)


directory2 <- paste(directory,"/Integration",sep="")

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
source("PSN_utils.R");
source("Diablo_utils.R");


# # Collect hyperparameters for the analysis (variables and values used defined)----
# int_method = "NEMO" # integration method (can be "SNF" or "NEMO")
# 
# K.nemo = NA # number of neighbors in KNN for NEMO (can also be integer or a vector)
# sigma = 0.5 # variance for affinityMatrix
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
# 

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
      directory2 <- paste(directory,"/Integration/INPUT/", "Metabolomics_", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
      serum_metabolomics <- vroom(paste(directory2,"/Metabolomics_",COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
      directory2 <- paste(directory,"/Metabolomics/OUTPUT/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
      serum_annotation <- vroom( paste(directory2,"/",COMMAND$LABEL[i],"_",args[1],"_vs_",args[2],"_results.tsv", sep = ""), delim = "\t")
      INPUTX<-Metabolomics_processing(serum_annotation,serum_metabolomics)
      assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
      myList <- list.append(myList,INPUTX[[1]])
      names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])    
    }
    
    if(COMMAND$DATA_TYPE[i] == "Transcriptomics"){
      
      print(args[1])
      directory2 <- paste(directory,"/Integration/INPUT/", COMMAND$LABEL[i],"_",args[1],"_vs_", args[2], sep ="")
      Wholeblood_RNAseq <-  vroom(paste(directory2, "/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], "_normalized_vst_variance.tsv",sep = ""), delim = "\t") #read normalization only
      Wholeblood_metadata <-  vroom(paste(directory2, "/","/Metadata_",COMMAND$LABEL[i], "_", args[1],".tsv",sep = ""), delim = "\t")
      INPUTX<-Transcriptomics_processing(Wholeblood_metadata,Wholeblood_RNAseq, Max_features_SNF)
      assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
      myList <- list.append(myList,INPUTX[[1]])
      names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])    
      
    }
    
    if(COMMAND$DATA_TYPE[i] == "Methylomics"){
      
      
      directory2 <- paste(directory,"/Integration/INPUT/", "Methylome_",COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="") 
      Methylome_WB <-  vroom(paste(directory2, "/", COMMAND$LABEL[i], "_matrix_MOFA.tsv",sep = ""), delim = "\t") #read normalization only
      Methylome_metadata <-  vroom(paste(directory2, "/", COMMAND$LABEL[i],"_metadata_MOFA.tsv",sep = "") ,delim = "\t")
      directory2 <- paste(directory,"/Methylomics/OUTPUT/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
      Methylome_annotation <- vroom(paste(directory2, "/", "DMP_", COMMAND$LABEL[i], "_Methylome_", args[1] ,"_vs_", args[2],".tsv",sep = ""), delim = "\t", col_names = TRUE)
      INPUTX<-Methylomics_processing(Methylome_annotation,Methylome_WB,Methylome_metadata, Max_features_SNF)
      assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
      myList <- list.append(myList,INPUTX[[1]])
      names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])    
    }
    
    if(COMMAND$DATA_TYPE[i] == "Undefined"){        
      directory2 <- paste(directory,"/Integration/INPUT/", "Undefined_", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
      samples_undefined <- vroom(paste(directory2,"/Undefined_",COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
      INPUTX<-Undefined_processing(samples_undefined)
      assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
      myList <- list.append(myList,INPUTX[[1]])
      names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])
      
    }        
  }
}

# Create directory to save all results ---
dir.create(path = paste(directory,"/Integration/OUTPUT/",int_method, "_", args[1] ,"_vs_", args[2],sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")

directory2 <- paste(directory,"/Integration/OUTPUT/",int_method, "_", args[1] ,"_vs_", args[2],sep="") 
setwd(directory2)

#Provide names to the list
names(myList) <- names_X

apply_scaling_SNF = as.logical(TRUE) # JG: it could be set as "TRUE"/"FALSE"
only_common_samples = as.logical(TRUE) # JG: it could be set as "TRUE"/"FALSE"
# Prepare data for NEMO----
# Feature selection used only to test the code, not necessary in final script
data <- snf_nemo.preprocess(myList, METADATA, 
                            fsel=TRUE, Max_features = Max_features_SNF, 
                            apply_scaling=apply_scaling_SNF, 
                            int_method = int_method, 
                            only_common_samples=only_common_samples)


# Set ground truth clustering (if available) ----
# This should be selected by the user

# get unified metadata dataframe
metadata_unified <- data$metadata

if(!is.na(gt.clust_name) && gt.clust_name != "X"){
    gt.clust <- as.factor(metadata_unified[[gt.clust_name]]) # factor (ground truth clustering)
    names(gt.clust) <- metadata_unified$ID
    
} else {
    gt.clust <- NULL # no ground truth clustering available
}

if (!gt.clust_name %in% colnames(data$metadata)) {
    stop("Ground-truth clustering variable not found in metadata: ", gt.clust_name)
}

# Check percentage of missing samples
if(!only_common_samples){
    
    missing_omics_summary <- warn_missing_omics_samples(
        data_list = data$data,
        metadata = metadata_unified,
        threshold = 0.20
    )
    
    write.csv(
        missing_omics_summary,
        file = "missing_samples_per_omic_NEMO.csv",
        row.names = FALSE
    )
    
}

# Define K for NEMO---
# if K.nemo is NA, it is computed as number of samples / NUM.NEIGHBORS.RATIO for
# each omic
if (any(is.na(K.nemo))) {
    K.nemo_vec = as.numeric(lapply(1:length(data$data), function(i) round(nrow(data$data[[i]]) / NUM.NEIGHBORS.RATIO)))
    warning("K.nemo is NA, it will be set to (number of samples / 6) for each omic.")
} else if (length(K.nemo) == 1) {
    K.nemo_vec = rep(K.nemo, length(data$data))
} else{
    #K.nemo_vec = K.nemo
    stop("K.nemo should be either NA or a single value.")
}


# Compute similarity matrix for each modality----
sim_mat <- lapply(seq_along(data$data), function(i) {
    scaled.exp.euclidean(data$data[[i]], K = K.nemo_vec[i], sigma = 0.5)
    })


names(sim_mat) <- names(data$data)

# save similarity matrices
dir.create("similarity_matrices", showWarnings = FALSE)
for(i in 1:length(sim_mat)){
    write.csv(sim_mat[[i]], file=file.path("similarity_matrices", 
              paste0(names(sim_mat)[i], "_similarity_matrix.csv")), 
              row.names=T)
}

# Check order of samples in metadata and similarity matrices
# check_names(sim_mat, data$metadata)

# Integrate using NEMO----
data_t <- lapply(data$data, t)
W_int <- nemo.affinity.graph_fix(data_t, k=K.nemo_vec)

# Reorder integrated similarity matrix and input data to match metadata
W_int <- W_int[match(data$metadata$ID, rownames(W_int)), match(data$metadata$ID, colnames(W_int))]

# Save integrated similarity matrix
write.csv(W_int, file=file.path("similarity_matrices", 
                                "integrated_similarity_matrix.csv"), row.names=T)

# Reorder gt.clust to match the integrated similarity matrix
if(!is.null(gt.clust)){
    gt.clust <- gt.clust[rownames(W_int)]
}

# Cluster optimization using eigengaps approach----

# Ground truth clustering (should be selected by user or set to NULL)
nc_estim <- estimate.nc(W_int, nc=nc, gt.clust=gt.clust, int_method=int_method)

# save optimal number of clusters
dir.create("clustering_metrics", showWarnings = FALSE)
write.csv(nc_estim, file="clustering_metrics/optimal_number_of_clusters.csv", 
          row.names=F)

# Spectral clustering----

# Apply clustering to all the number of clusters in nc_range
# Note that NEMO uses the spectralClustering implemented by SNFtool,
# see github (https://github.com/Shamir-Lab/NEMO/blob/451052feb98cd58dcc7e3db2171217a95066fb69/NEMO/R/NEMO.R#L54)
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
if(is.na(K.nemo)){
    K_W_int <- round(ncol(W_int) / NUM.NEIGHBORS.RATIO)
} else if(length(K.nemo) == 1){
    K_W_int <- K.nemo
}

# K for each modality and integrated matrix
K_sparse <- c(K.nemo_vec, K_W_int) 


make_graph(c(sim_mat, list(W_int=W_int)), gt.clust=gt.clust, 
           pred.clust=pred.clust, sparse_method="KNN", K=K_sparse, path="graphs")

# # Compute external quality indices for clustering----
# # External quality indices can be computed only if a ground truth clustering is provided!!!
types <- c("ARI", "AMI", "NVI", "NMI")
if(!is.null(gt.clust)){
    y = as.numeric(gt.clust)
    ext.val.idx <- apply(clustering, 2,
                         function(x) {
                             res <- as.data.frame(aricode::clustComp(x, y)[types]);
                             #res$"V-measure" <- sabre::vmeasure(x, y)$v_measure
                             rownames(res) <- NULL;
                             res})
    for(kn in names(ext.val.idx)){
        ext.val.idx[[kn]]$nc <- kn
    }
    ext.val.idx <- do.call(rbind, ext.val.idx)
    rownames(ext.val.idx) <- NULL
}

# # Visualize and save external quality indices by number of clusters----
# # External quality indices can be computed only if a ground truth clustering is provided!!!
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
  METADATA_alluvial <-  as.data.frame(data$metadata[, valid_enrich_vars, drop=F])
  rownames(METADATA_alluvial) <- data$metadata$ID
  
  plot.alluvial(clustering[,as.character(nc_estim$nc_estim), drop=F], gt.clust, 
                METADATA_alluvial, 
                save_path="./clustering_metrics")

}

# Which are the most important features for integration?----
feat_imp <- compute_feature_importance(data$data, pred.clust, K=K.nemo_vec, sigma=0.5)

if (!is.null(feat_imp)) {
  dir.create("feature_importance", showWarnings = FALSE)
  write.csv(feat_imp, file="feature_importance/feature_importance.csv", row.names=F)
  
  # Heatmap of top_feat important features
  plot_imp_heatmap(data$data, feat_imp, pred.clust=pred.clust, gt.clust=gt.clust, 
                   top_feat=top_feat, save_path="./feature_importance")
}

# Compute variables enrichment and log-rank test----

# add fake survival data (just to test the function)
# set.seed(123)
# metadata <- metadata_unified # identical in different modalities
# metadata$OS.time <- sample(1:100, nrow(metadata), replace=T)
# metadata$OS.event <- sample(0:1, nrow(metadata), replace=T)
# metadata$DFS.time <- sample(1:100, nrow(metadata), replace=T)
# metadata$DFS.event <- sample(0:1, nrow(metadata), replace=T)
# enrich_vars <- c("GENDER", "AGE")
# surv_vars <- c("OS.time", "OS.event", "DFS.time", "DFS.event")


dir.create("enrichment_survival_analysis", showWarnings = FALSE)

valid_surv_vars <- surv_vars[surv_vars != "X"]

if (length(valid_surv_vars) > 0 || length(valid_enrich_vars) > 0) {
  
  #quiet_enrich_surv_analysis <- quietly(enrich_surv_analysis) # Catch warnings
  enrich_surv_res <- enrich_surv_analysis(clustering, metadata=data$metadata, 
                                          valid_enrich_vars, valid_surv_vars, 
                                          file_path="enrichment_survival_analysis")
  
  enrich_surv_res_save <- rbind(enrich_surv_res$result$enrich_res, 
                                enrich_surv_res$result$surv_res)
  write.csv(enrich_surv_res_save, 
            file="enrichment_survival_analysis/enrichment_survival_results.csv", 
            row.names=T)
  
  write.csv(enrich_surv_res$warning_log, 
            file="enrichment_survival_analysis/enrichment_survival_results_warnings.csv", 
            row.names=T)
}

