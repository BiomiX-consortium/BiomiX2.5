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


# # #MANUAL INPUT
args = as.list(c("BLymphocytes","SJS"))
args[2] <-"unmutated"
args[1] <-"mutated"
args[3] <- "C:/Users/crist/Desktop/BiomiX2.5"

directory <-args[3]
# 

int_method <- "SNF"

MART <- vroom(paste(directory,"/Integration/x_BiomiX_DATABASE/mart_export_37.txt",sep=""), delim = ",")
myList <- list()

COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")
COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
Max_features <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[3])
Max_features_SNF <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[3])


K.snf <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[1]) # number of neighbors in KNN
sigma <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[2]) # variance for affinityMatrix
t <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS[3]) # number of iterations for SNF
nc <- c(2:as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[1])) # Max number of cluster to test
top_feat <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS[2]) #Number variable in the heatmap to visualize
# Variable of interest for enrichment and survival analysis
enrich_vars <- strsplit(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[1]), "/")[[1]]
surv_vars <- strsplit(as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[2]), "/")[[1]]
DIR_METADATA <- readLines(paste(directory,"directory.txt",sep="/"))

# Ground truth clustering name (if available, otherwise NULL)
gt.clust_name <- as.character(COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES[3])



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

# Read data----
# For now, I use example data in "Integration_data" folder. 
# In the future, data normalized by unimodal biomix pipelines will be used

# transcriptomics + methylomics
# input_path = "C:/Users/crist/Downloads/Integration_data/Integration_data/CLL"; 

# transcriptomics + metabolomics (undefined)
#input_path = "../../Integration_data/Tubercolosis" 


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
    INPUTX<-Metabolomics_processing(serum_annotation,serum_metabolomics,COMMAND$LABEL[i])
    assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
    myList <- list.append(myList,INPUTX[[1]])
    names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])    
  }
  
  if(COMMAND$DATA_TYPE[i] == "Transcriptomics"){
    
    print(args[1])
    directory2 <- paste(directory,"/Integration/INPUT/", COMMAND$LABEL[i],"_",args[1],"_vs_", args[2], sep ="")
    Wholeblood_RNAseq <-  vroom(paste(directory2, "/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], "_normalized_vst_variance.tsv",sep = ""), delim = "\t") #read normalization only
    Wholeblood_metadata <-  vroom(paste(directory2, "/","/Metadata_",COMMAND$LABEL[i], "_", args[1],".tsv",sep = ""), delim = "\t")
    INPUTX<-Transcriptomics_processing(Wholeblood_metadata,Wholeblood_RNAseq,COMMAND$LABEL[i])
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
    INPUTX<-Methylomics_processing(Methylome_annotation,Methylome_WB,Methylome_metadata,COMMAND$LABEL[i])
    assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
    myList <- list.append(myList,INPUTX[[1]])
    names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])    
  }
  
  if(COMMAND$DATA_TYPE[i] == "Undefined"){        
    directory2 <- paste(directory,"/Integration/INPUT/", "Undefined_", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
    samples_undefined <- vroom(paste(directory2,"/Undefined_",COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
    INPUTX<-Undefined_processing(samples_undefined,COMMAND$LABEL[i])
    assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
    myList <- list.append(myList,INPUTX[[1]])
    names_X <- list.append(names_X,COMMAND$DATA_TYPE[i])
    
  }        
}
}

# Create directory to save all results ---
# timestamp <- format(Sys.time(), "%d-%m-%Y_%H-%M-%S")
# save_path <- file.path(basename(input_path), paste(int_method, timestamp, sep="_"))
# dir.create(save_path, showWarnings = FALSE, recursive=TRUE)
# setwd(save_path); # work inside data folder

dir.create(path = paste(directory,"/Integration/OUTPUT/",int_method, "_", args[1] ,"_vs_", args[2],sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")

directory2 <- paste(directory,"/Integration/OUTPUT/",int_method, "_", args[1] ,"_vs_", args[2],sep="") 
setwd(directory2)

#Provide names to the list
names(myList) <- names_X


# Prepare data for SNF----
# Feature selection used only to test the code, not necessary in final script
data <- snf.preprocess(myList, METADATA,
                       fsel=TRUE, Max_features = Max_features_SNF)


#X Jessica: The feature selection and the matrix transformation are already occurring before
# data <- snf.preprocess(myList, METADATA,
#                        fsel=TRUE, Max_features = 200)

# Set ground truth clustering (if available) ----
# This should be selected by the user
if(!is.null(gt.clust_name)){
    # metadata are the same for all omics
    gt.clust <- as.factor(data$metadata[[gt.clust_name]]) # factor (ground truth clustering)
    names(gt.clust) <- data$metadata$ID
} else {
    gt.clust <- NULL # no ground truth clustering available
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
check_names(sim_mat, data$metadata, int_method = int_method)

print("Dimension check: The matrix must be of the same size")
lapply(sim_mat, dim)  # Check dimensions

# Integrate using SNF or NEMO----
W_int <- SNFtool::SNF(sim_mat, K=K.snf, t=t)


# Save integrated similarity matrix
write.csv(W_int, file=file.path("similarity_matrices", 
                                "integrated_similarity_matrix.csv"), row.names=T)


# Cluster optimization using eigengaps approach----

# Ground truth clustering (should be selected by user or set to NULL)
nc_estim <- estimate.nc(W_int, nc=nc, gt.clust=gt.clust, int_method=int_method)

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
# External quality indices can be computed only if a ground truth clustering is provided!!!
# types <- c("ARI", "AMI", "NVI", "NMI")
# if(!is.null(gt.clust)){
#     y = as.numeric(gt.clust)
#     ext.val.idx <- apply(clustering, 2, 
#                          function(x) {
#                              res <- as.data.frame(aricode::clustComp(x, y)[types]);
#                              res$"V-measure" <- sabre::vmeasure(x, y)$v_measure
#                              rownames(res) <- NULL;
#                              res})
#     for(kn in names(ext.val.idx)){
#         ext.val.idx[[kn]]$nc <- kn
#     }
#     ext.val.idx <- do.call(rbind, ext.val.idx)
#     rownames(ext.val.idx) <- NULL
# }
# 
# 
# # Visualize and save external quality indices by number of clusters----
# # External quality indices can be computed only if a ground truth clustering is provided!!!
# if(!is.null(gt.clust)){
#     plot_ext.val.idx(ext.val.idx, nc_estim$nc_range)
# }

print(enrich_vars)
METADATA_alluvial <- data$metadata[, enrich_vars, drop=F]
rownames(METADATA_alluvial) <- data$metadata$ID

# Make Alluvial plot to compare with ground truth clustering and variables----
# For now, I compare the optimal cluster to gt.clust (if present) and variables of interest (if present).
plot.alluvial(clustering[,as.character(nc_estim$nc_estim), drop=F], gt.clust, 
              METADATA_alluvial, 
              save_path="./clustering_metrics")

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
set.seed(123)
metadata <- METADATA # identical in different modalities
metadata$OS.time <- sample(1:100, nrow(metadata), replace=T)
metadata$OS.event <- sample(0:1, nrow(metadata), replace=T)
metadata$DFS.time <- sample(1:100, nrow(metadata), replace=T)
metadata$DFS.event <- sample(0:1, nrow(metadata), replace=T)

dir.create("enrichment_survival_analysis", showWarnings = FALSE)

if (surv_vars != "X") {

quiet_enrich_surv_analysis <- quietly(enrich_surv_analysis) # Catch warnings
enrich_surv_res <- quiet_enrich_surv_analysis(clustering, metadata=metadata, 
                                        enrich_vars, surv_vars, 
                                        file_path="enrichment_survival_analysis")

enrich_surv_res_save <- rbind(enrich_surv_res$result$enrich_res, 
                              enrich_surv_res$result$surv_res)
write.csv(enrich_surv_res_save, 
          file="enrichment_survival_analysis/enrichment_survival_results.csv", 
          row.names=T)

write.csv(enrich_surv_res$warnings, 
          file="enrichment_survival_analysis/enrichment_survival_results_warnings.csv", 
          row.names=T)
}

