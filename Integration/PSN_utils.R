# Functions used in main script "SNF_int.R"
# Author: Jessica Gliozzo


#' Read data and metadata from TSV files
#'
#' @param input_path string. Path to directory containing data and metadata. 
#' Each modality is in a separate folder, named by the modality type, and it 
#' contains two files (data and metadata). 
#' Note that "data" and "metadata" have to be part of the filename.
#'
#' @return A list with two elements:
#' - Data: a list with a dataframe for each modality.
#' - Metadata: a list with a dataframe for each modality.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
read_data <- function(input_path){
    
    folders <- list.dirs(input_path, recursive = FALSE)
    
    data = list()
    for(folder in folders){
        
        files <- list.files(folder, full.names = TRUE)
        for(file in files){
            
            # Read tsv file of data and metadata
            path_split <- tolower(strsplit(file, "/")[[1]])
            omic_name <- path_split[[length(path_split)-1]]
            data_type <- ifelse(grepl("metadata", path_split[[length(path_split)]]), "metadata", "data")
            
            data[[data_type]][[omic_name]] <- read.table(file = file, sep = '\t', header = TRUE)
        }
    }
    
    return(data)
}


#' Prepare data and metadata for SNF integration
#' 
#' @description
#' This function performs the following pre-processing steps:
#' - Remove features having missing values
#' - Feature selection by variance (optional)
#' - Z-score standardization
#' - Retain the subset of samples having all omics in data and metadata
#' 
#'
#' @param data List. List of dataframes (features x samples), one for each modality.
#' @param metadata List. List of dataframes (samples x features), one for each modality.
#' @param col.remove vector. Vector containing the names of columns in data and/or 
#' metadata to remove. 
#' @param fsel logical. Perform feature selection by variance.
#' @param Max_features numeric. Number of features to be retained in the 
#' filter by variance (used only if fsel=TRUE).
#'
#' @return A list with two elements:
#' - Data: a list with a dataframe for each modality (samples x features).
#' - Metadata: a list with a dataframe for each modality (samples x features).
#' 
#' @author Jessica Gliozzo
#' 
#' @export
snf.preprocess <- function(data, metadata,
                           fsel=FALSE, Max_features = 200){
    
    for(omic in names(data)){
        
        if(omic == "Methylomics"){
            feature_names <- rownames(data[[omic]])
        } else if (omic == "Transcriptomics"){
            feature_names <- rownames(data[[omic]])
        } else if (omic == "Metabolomics" | omic == "Undefined"){
          feature_names <- rownames(data[[omic]])
        }
        
        else {
            stop("Unknown omic type!")
        }
        
        # Retain features having no missing values
        keep <- complete.cases(data[[omic]])
        data[[omic]] <- data[[omic]][keep, ]
        feature_names <- feature_names[keep]
        
        # Convert to numeric
        data[[omic]] <- apply(data[[omic]], 2, as.numeric)
        
        # Set features names
        rownames(data[[omic]]) <- feature_names
        
        # Feature selection by variance and transpose data (sample x features)
        if(fsel){
            var_filter <- function(mat, Max_features){
                # Stop if the number of features is less than the number of features to retain
                if(ncol(mat) < Max_features){
                    stop(paste0(omic, ": The number of features is less than the number of features to retain! Did you select more SNF integration variables than the general integration input features?"))
                }
                
                variances <- apply(mat, 2, var)
                sorted <- sort(variances, decreasing=TRUE, index.return=TRUE)$ix[1:Max_features]
                
                return(mat[,sorted])
            }
            
            data[[omic]] <- var_filter(data[[omic]], Max_features)
        } else {
            print("No SNF Filtering")
        }
        
        # Add warning if the number of features is greater than 5000
        if(ncol(data[[omic]]) > 5000){
            warning(paste0(omic, ": The number of features is greater than 5000! Consider feature selection."))
        }
        
      
        #Standardize data (already done in the single omics pipeline) *Cristian
        # data[[omic]] <- SNFtool::standardNormalization(data[[omic]])
        
    }
    
    # get common samples among all omics
    common_samples <- Reduce(intersect, lapply(data, rownames))
    metadata <- metadata[which(metadata$ID %in% common_samples),]

        # Collect common set of samples in all omics since SNF does not natively
    # handle missing data
    for(omic in names(data)){
        data[[omic]] <- data[[omic]][common_samples, ]
        
    }
    
    return(list(data=data, metadata=metadata))
}


#' Prepare data and metadata for NEMO integration
#' 
#' @description
#' This function performs the following pre-processing steps:
#' - Remove features having missing values
#' - Feature selection by variance (optional)
#' - Z-score standardization
#' 
#' NEMO is able to integrate data from different omics where some samples are missing.
#'
#' @param data List. List of dataframes (features x samples), one for each modality.
#' @param metadata List. List of dataframes (samples x features), one for each modality.
#' @param col.remove vector. Vector containing the names of columns in data and/or 
#' metadata to remove. 
#' @param fsel logical. Perform feature selection by variance.
#' @param Max_features numeric. Number of features to be retained in the 
#' filter by variance (used only if fsel=TRUE).
#'
#' @return A list with two elements:
#' - Data: a list with a dataframe for each modality (features x samples).
#' - Metadata: a list with a dataframe for each modality (samples x features).
#' 
#' @author Jessica Gliozzo
#' 
#' @export
nemo.preprocess <- function(data, metadata,
                           fsel=FALSE, Max_features = 200){
  
  for(omic in names(data)){
    
    if(omic == "methylomics"){
      feature_names <- data[[omic]]$ID
    } else if (omic == "transcriptomics"){
      feature_names <- data[[omic]]$Gene.name
    } else if (omic == "metabolomics" | omic == "undefined"){
      rownames(data[[omic]]) <- data[[omic]]$ID
    }
    
    else {
      stop("Unknown omic type!")
    }
    
    # Remove unwanted columns
    data[[omic]] <- as.matrix(data[[omic]][, !names(data[[omic]]) %in% col.remove])
    
    # Transpose data to keep the same preprocessing steps
    if(omic == "metabolomics" | omic == "undefined"){
      data[[omic]] <- t(data[[omic]])
      feature_names <- rownames(data[[omic]])
    }
    
    # Retain features having no missing values
    keep <- complete.cases(data[[omic]])
    data[[omic]] <- data[[omic]][keep, ]
    feature_names <- feature_names[keep]
    
    # Convert to numeric
    data[[omic]] <- apply(data[[omic]], 2, as.numeric)
    
    # Set features names
    rownames(data[[omic]]) <- feature_names
    
    # Feature selection by variance
    if(fsel){
      var_filter <- function(mat, Max_features){
        # Stop if the number of features is less than the number of features to retain
        if(nrow(mat) < Max_features){
          stop(paste0(omic, ": The number of features is less than the number of features to retain!"))
        }
        
        variances <- apply(mat, 1, var)
        sorted <- sort(variances, decreasing=TRUE, index.return=TRUE)$ix[1:Max_features]
        
        return(mat[sorted, ])
      }
      
      data[[omic]] <- var_filter(data[[omic]], Max_features)
    }
    
    # Add warning if the number of features is greater than 5000
    if(nrow(data[[omic]]) > 5000){
      warning(paste0(omic, ": The number of features is greater than 5000! Consider feature selection."))
    }
    
    
    #Standardize data
    data[[omic]] <- t(SNFtool::standardNormalization(t(data[[omic]])))
    
    # Convert to data.frame
    data[[omic]] <- as.data.frame(data[[omic]])
  }
  
  # Fix metadata
  for(omic in names(data)){
    
    if(!is.null(metadata[[omic]])){
      # subset metadata
      rownames(metadata[[omic]]) <- metadata[[omic]]$ID
      
      # Fix column data type in metadata
      metadata[[omic]] <- metadata[[omic]] %>% 
        dplyr::mutate(across(!ID & where(is.character), as.factor))
    }
  }
  
  # Check that data and metadata have the same samples in the same order
  for(omic in names(data)){
    if(!is.null(metadata[[omic]])){
      if(!all(colnames(data[[omic]]) == rownames(metadata[[omic]]))){
        stop(paste0("Data and metadata for ", omic, " do not have the same samples in the same order!"))
      }
    } 
  }
  
  return(list(data=data, metadata=metadata))
}


#' Compute the scaled exponential euclidean kernel
#'
#' @param M dataframe/matrix.
#' @param K integer. Number of nearest neighbours.
#' @param sigma integer. Variance for local model.
#'
#' @return Affinity matrix (samples x samples).
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#'
scaled.exp.euclidean <- function(M, K=20, sigma=0.5){
    
    dist <- (SNFtool::dist2(as.matrix(M), as.matrix(M)))^(1/2);
    sim <- SNFtool::affinityMatrix(dist, K=K, sigma=sigma);
    
    return(sim)
}


#' Plot heatmaps of affinity matrices
#'
#' @param mat_list list. Named list of affinity matrices.
#' @param gt.clust factor. Factor containing the ground truth clustering (def. NULL).
#' @param pred.clust vector. Named vector of predicted clustering.
#' @param path string. Path to save results.
#' @param norm string. Type of normalization to apply on the affinity matrix.
#' Possible options are: "row" (normalize by sum of rows), "max" (max normalization), 
#' "none" (no normalization).
#'
#' @return Plot an heatmap for each affinity matrix.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#'
make_heatmap <- function(mat_list, gt.clust=NULL, pred.clust, 
                         path=".", norm="row"){
    
    # Get df with ground truth clustering
    if(!is.null(gt.clust)){
        ind <- sort(as.vector(gt.clust), index.return=TRUE)
        ind <- ind$ix
        
        gt.clust.df <- data.frame("Clusters"=gt.clust[ind])
        rownames(gt.clust.df) <- names(gt.clust)[ind]
        
        idx <- match(rownames(gt.clust.df), names(pred.clust))
        gt.clust.df[["Pred. Clusters"]] <- as.factor(pred.clust[idx])
    }    
    
    for(omic in names(mat_list)){
        mat <- mat_list[[omic]]
        
        #Normalize matrix
        if(norm == "row"){
            mat <- W_normalize(mat)
        } else if(norm == "max"){
            diag(mat) <- 0
            mat <- mat/max(mat) #0-1 normalization
        } else if(norm == "none"){
            mat <- mat
        } else {
            stop("Unknown normalization method!")
        }
        
        if(omic == "W_int"){omic <- "integrated_data"}
        
        # Make plot
        if(!is.null(gt.clust)){
          
            # subset matrix and ground truth clustering by common samples
            common_samples <- intersect(rownames(gt.clust.df), rownames(mat))
            mat_sub <- mat[common_samples, common_samples]
            gt.clust.df_sub <- gt.clust.df[common_samples, , drop=FALSE]
          
            # plot
            png(file.path(path, paste0("Similarity_heatmap_",omic, ".png")), 
                width=4, height=4, units = "in", res=300)
            pheatmap::pheatmap(mat_sub, cluster_rows = F, 
                               cluster_cols = F, annotation_col = gt.clust.df_sub, 
                               show_colnames = F, show_rownames = F, 
                               fontsize = 6)
            dev.off()
            
        } else {
            
            # order by predicted clustering
            ind <- sort(as.vector(pred.clust), index.return=TRUE)
            ind <- ind$ix
            pred.clust.df <- data.frame(pred.clust=as.factor(pred.clust[ind]))
            rownames(pred.clust.df) <- names(pred.clust)[ind]
            
            # subset matrix and predicted clustering by common samples
            common_samples <- intersect(rownames(pred.clust.df), rownames(mat))
            mat_sub <- mat[common_samples, common_samples]
            pred.clust.df_sub <- pred.clust.df[common_samples, , drop=FALSE]
            
            # plot
            png(file.path(path, paste0("Similarity_heatmap_",omic, ".png")),
                width=4, height=4, units = "in", res=300)
            pheatmap::pheatmap(mat_sub, cluster_rows = F, 
                               cluster_cols = F, annotation_col = pred.clust.df_sub, 
                               show_colnames = F, show_rownames = F, fontsize = 6)
            dev.off()
        }
        
    } 
    
}


#' Normalize affinity matrix
#' 
#' @description
#' Normalize matrix by sum of its rows.
#' 
#'
#' @param W matrix.
#'
#' @return Normalized matrix.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#' 
W_normalize <- function(W){
    
    diag(W) <- 0
    W <- W / rowSums(W)
    W <- W + t(W)
    
    return(W)
}


#' Apply k-Nearest Neighbour sparsification
#' 
#' @description
#' Retain only the first k neighbous of each sample.
#' 
#'
#' @param sim.mat matrix. Affinity matrix.
#' @param k integer. Number of neighbours to retain.
#'
#' @return A sparsified symmetric affinity matrix.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#'
knn_sparsification <- function(sim.mat, k){
    
    # Do nor consider self-loops
    diag(sim.mat) <- 0
    
    # Compute knn of each row
    compute.row.knn <- function(sim.row, k=k) {
        knn.row = sim.row
        thresh = sort(sim.row, decreasing = T)[k]
        knn.row[sim.row < thresh] = 0
        row.sum = sum(knn.row)
        knn.row[sim.row >= thresh] = knn.row[sim.row >= thresh] / row.sum
        
        return(knn.row)
    }
    
    
    mat.knn <- apply(sim.mat, 1, compute.row.knn, k=k) # apply to each row and saved in columns
    
    # Make symmetric
    mat.knn <- mat.knn + t(mat.knn)
    #mat.knn <-  pmax(mat.knn, t(mat.knn))
    
    return(mat.knn)
}


# Define a custom triangle shape
triangle_shape <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    vertex.size  <- params("vertex", "size")
    
    for(i in 1:nrow(coords)){
        x <- coords[i, 1]
        y <- coords[i, 2]
        # Use vertex.size to scale the triangle (adjust divisor to suit your plotting scale)
        radius <- vertex.size[i] / 200  
        angles <- c(pi/2, pi/2 + 2*pi/3, pi/2 + 4*pi/3)
        xs <- x + radius * cos(angles)
        ys <- y + radius * sin(angles)
        polygon(xs, ys, col = vertex.color[i], border = "black")
    }
}
add_shape("triangle", clip = igraph.shape.noclip, plot = triangle_shape)

# Define a custom star shape (5-point star)
star_shape <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    vertex.size  <- params("vertex", "size")
    
    for(i in 1:nrow(coords)){
        x <- coords[i, 1]
        y <- coords[i, 2]
        radius <- vertex.size[i] / 200
        
        # Outer angles for the 5 points
        outer_angles <- seq(pi/2, 2*pi + pi/2, length.out = 6)[-6]
        # Inner angles: offset by pi/5 (36°) relative to the outer points
        inner_angles <- outer_angles + pi/5
        
        # Interleave outer and inner angles
        all_angles <- c(rbind(outer_angles, inner_angles))
        # Alternate radii: outer radius for the outer points, inner radius (half of outer) for the inner points
        all_radii <- rep(c(radius, radius/2), length.out = length(all_angles))
        
        xs <- x + all_radii * cos(all_angles)
        ys <- y + all_radii * sin(all_angles)
        
        polygon(xs, ys, col = vertex.color[i], border = "black")
    }
}
add_shape("star", clip = igraph.shape.noclip, plot = star_shape)

# Define a custom diamond shape
diamond_shape <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    vertex.size  <- params("vertex", "size")
    
    for(i in 1:nrow(coords)){
        x <- coords[i, 1]
        y <- coords[i, 2]
        radius <- vertex.size[i] / 200
        xs <- c(x, x + radius, x, x - radius)
        ys <- c(y + radius, y, y - radius, y)
        polygon(xs, ys, col = vertex.color[i], border = "black")
    }
}
add_shape("diamond", clip = igraph.shape.noclip, plot = diamond_shape)

# Define a custom hexagon shape
hexagon_shape <- function(coords, v = NULL, params) {
    vertex.color <- params("vertex", "color")
    vertex.size  <- params("vertex", "size")
    
    for(i in 1:nrow(coords)) {
        x <- coords[i, 1]
        y <- coords[i, 2]
        radius <- vertex.size[i] / 200
        # Generate 6 angles for the hexagon (0 to 2*pi, excluding the duplicate final point)
        angles <- seq(0, 2*pi, length.out = 7)[-7]
        xs <- x + radius * cos(angles)
        ys <- y + radius * sin(angles)
        polygon(xs, ys, col = vertex.color[i], border = "black")
    }
}
add_shape("hexagon", clip = igraph.shape.noclip, plot = hexagon_shape)


#' Plot network from affinity matrix
#'
#' @param mat_list list. Named list of affinity matrices.
#' @param gt.clust factor. Factor containing the ground truth clustering (def. NULL).
#' @param pred.clust vector. Named vector of predicted clustering.
#' @param path string. Path to save results.
#' @param sparse_method string. Method to sparsity the graph for visualization.
#' Options are "KNN" (use KNN sparsification), "quantile" (cut edges having weight
#' below second quartile - median) and "none" (no sparsification).
#' @param K vector. vector with the number of neighbours to retain (used only if 
#' sparse_method=="KNN") in each affinity matrix.
#'
#' @return Return a static (.png) and interactive (.html) plot for each affinity matrix.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#'
make_graph <- function(mat_list, gt.clust=NULL, pred.clust, path=".", 
                       sparse_method="KNN", K){
    
    pal <- c("lightgreen", "plum", "darkorange", "gold", "lightblue", "forestgreen", 
             "firebrick", "blue")
    
    count=1
    for(omic in names(mat_list)){
        
        if(omic == "W_int"){
            omic_name <- "integrated data"
        } else {
            omic_name <- omic
        }
        
        main_title <- tools::toTitleCase(paste("similarity graph for", omic_name))
        
        mat <- mat_list[[omic]]
        
        # subset matrix by common samples in pred.clust and gt.clust (if present)
        if(!is.null(gt.clust)){
            common_samples <- intersect(rownames(mat), names(gt.clust))
            mat <- mat[common_samples, common_samples]
            gt.clust <- gt.clust[common_samples]
            pred.clust <- pred.clust[common_samples]
        } else {
            common_samples <- intersect(rownames(mat), names(pred.clust))
            mat <- mat[common_samples, common_samples]
            pred.clust <- pred.clust[common_samples]
        }
        
        #Normalize matrix
        if(sparse_method == "KNN"){
            # sparsify graph using KNN
            mat <- knn_sparsification(mat, k=K[count])
            message(paste0("Sparsification method for ", omic,": KNN with k=", K[count]))
        } else if(sparse_method == "quantile"){
            # cut edges of similarity matrix based on quantiles (0.5%)
            q.value <- quantile(mat, probs = c(.25, .5, .75))[2]
            mat[mat < q.value] <- 0
        } else if(sparse_method == "none"){
            mat <- mat
        } else {
            stop("Unknown sparsification method!")
        }
        
        # Make plot
        if(!is.null(gt.clust)){
            # Color nodes by gt clustering and shape by predicted clustering
            
            g <- igraph::graph_from_adjacency_matrix(mat, mode="undirected", 
                                                     weighted=TRUE, diag=FALSE)
            
            
            possible.shapes <- c("circle", "square", "star", "triangle", "diamond", "hexagon")
            if(length(unique(pred.clust)) > length(possible.shapes)){
                warning("Number of predicted clusters is greater than the number of possible shapes! 
                        Some clusters will have the same shape.")
            }
            possible.shapes <- rep(possible.shapes, ceiling(length(unique(pred.clust))/length(possible.shapes)))
            shapes <- possible.shapes[1:length(unique(pred.clust))]
            V(g)$shape <- shapes[as.factor(pred.clust)]
            
            set.seed(123)
            if(length(unique(gt.clust)) > 8){
                warning("Number of ground truth clusters is greater than the number of colors! 
                        Some clusters will have the same color.")
            }
            
            colors <- rep(pal, ceiling(length(unique(gt.clust))/8))
            colors <- colors[1:length(unique(gt.clust))]
            V(g)$color <- colors[gt.clust]
            
            png(file.path(path, paste0("Similarity_graph_",omic, ".png")), 
                width=1000, height=700, res=150)
            par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
            plot(g, vertex.label.cex=0.6, vertex.label.color="black", 
                 main=main_title)
            
            legend("topright", 
                   legend = levels(gt.clust),
                   pch = 21,
                   pt.bg = colors,
                   title = "Clusters",
                   bty = "n", inset=c(-0.25,0))
            
            pch_vec <- c(1, 0, 8, 2, 18, 72)
            pch_vec <- rep(pch_vec, ceiling(length(unique(pred.clust))/length(pch_vec)))
            pch_vec <- pch_vec[1:length(unique(pred.clust))]
            legend("bottomright", 
                   legend = levels(as.factor(pred.clust)),
                   pch = pch_vec,
                   title = "Pred. Clusters",
                   bty = "n", inset=c(-0.25,0))
            dev.off()
            
            shapes <- ifelse(shapes == "circle", "dot", shapes)
            
            nodes_attr <- data.frame(color=rep(colors, each=length(shapes)), 
                                     shape=rep(shapes, times=length(unique(colors))), 
                                     label=paste0("Clust. (color): ", rep(levels(gt.clust), 
                                                                       each=length(unique(pred.clust))), 
                                                  ";\n Pred. Clust. (shape):", rep(levels(as.factor(pred.clust)), 
                                                                         times=length(levels(gt.clust)))
                                     ))
            
            
            # save interactive graph
            V(g)$sel <- paste(V(g)$color, V(g)$shape, sep = ",")
            ig.save <- visIgraph(g) %>%
                visOptions(highlightNearest = list(enabled = T, hover = T), 
                           nodesIdSelection = list(enabled=T, main="Select by name"), 
                           selectedBy = list(variable="sel", main="Select by color or shape", multiple=T)) %>%
                visLegend(position = "right", useGroups=F, addNodes = nodes_attr, 
                          main = main_title, width = 0.1)  %>%
                htmlwidgets::onRender("function(el, x) {
                                            el.style.height = '600px';
                                            el.style.width = '100%';
                                      }")
            saveWidget(ig.save, file = file.path(path, paste0("Similarity_graph_",omic, ".html")),selfcontained = FALSE)
            
            
        } else {
            g <- igraph::graph_from_adjacency_matrix(mat, mode="undirected", 
                                                     weighted=TRUE, diag=FALSE)
            
            
            set.seed(123)
            if(length(unique(pred.clust)) > 8){
                warning("Number of clusters is greater than the number of colors! 
                        Some clusters will have the same color.")
            }
            colors <- rep(pal, ceiling(length(unique(pred.clust))/8))
            colors <- colors[1:length(unique(pred.clust))]
            V(g)$color <- colors[pred.clust]
            
            png(file.path(path, paste0("Similarity_graph_",omic, ".png")), 
                width=1000, height=700, res=150)
            par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
            plot(g, vertex.label.cex=0.6, vertex.label.color="black", 
                 main=main_title)
            
            legend("right", 
                   legend = levels(as.factor(pred.clust)),
                   pch = 21,
                   pt.bg = colors,
                   title = "Pred. Clusters",
                   bty = "n", inset=c(-0.25,0))
            dev.off()
            
            nodes_attr <- data.frame(label=paste("Pred.:", sort(unique(pred.clust))), 
                                     color=colors, 
                                     shape= rep("dot", length(unique(pred.clust))))
            
            # save interactive graph
            ig.save <- visIgraph(g) %>%
                visLegend(position = "right", useGroups=F, addNodes = nodes_attr,
                          main = main_title, width = 0.1) %>% 
                visOptions(highlightNearest = list(enabled = T, hover = T), 
                          nodesIdSelection = list(enabled=T, main="Select by name"), 
                          selectedBy = list(variable="color", main="Select by color")) %>%
                htmlwidgets::onRender("function(el, x) {
                                            el.style.height = '600px';
                                            el.style.width = '100%';
                                      }")
            saveWidget(ig.save, file = file.path(path, paste0("Similarity_graph_",omic, ".html")))
            
        }
      
      count <- count + 1  
    } 
    
}


#' Estimate the optimal number of clusters
#'
#' @description
#' The eigengap method is used for estimation if int_method is "SNF". 
#' A modified version of the eigengap method is used if int_method is "NEMO".
#' 
#'
#' @param W matrix. Affinity matrix (samples x samples).
#' @param nc vector. Number of clusters to evaluate.
#' @param gt.clust factor. Factor containing the ground truth clustering 
#' (def. NULL).
#' @param int_method string. Integration method used for clustering. If can 
#' be "SNF" or "NEMO".
#'
#' @return List containing:
#' - The optimal number of clusters "nc_estim" estimated by eigengap method.
#' - The considered range of values.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
estimate.nc <- function(W, nc=NULL, gt.clust=NULL, int_method){
    
    # Check if int_method is valid
    if(!(int_method %in% c("SNF", "NEMO"))){
        stop("Unknown integration method! Use 'SNF' or 'NEMO'.")
    }
    
    if(is.null(gt.clust)){
        # Case without ground truth clustering
        if(is.null(nc)){
            max_nc <- min(floor(nrow(W) / 5), 10)
            nc <- 2:max_nc
            if(int_method == "SNF"){
              nc_estim <- SNFtool::estimateNumberOfClustersGivenGraph(W_int, NUMC=nc)$`Eigen-gap best`
            } else if(int_method == "NEMO"){
              nc_estim <- NEMO::nemo.num.clusters(W_int, NUMC=nc)
            }
            
            warning("The number of clusters is automatically estimated in range 2:", max_nc, ". Consider to set the number of clusters based on your sample size!")
            
        } else if(length(nc) == 1){
            nc_estim <- nc
            message("The number of clusters is set to ", nc, "!")
            
        } else if(length(nc) > 1){
            if(int_method == "SNF"){
              nc_estim <- SNFtool::estimateNumberOfClustersGivenGraph(W_int, NUMC=nc)$`Eigen-gap best`
            } else if(int_method == "NEMO"){
              nc_estim <- NEMO::nemo.num.clusters(W_int, NUMC=nc)
            }
            warning("The number of clusters is automatically estimated in range: ", min(nc),":", max(nc))
        }
        
    } else {
        # Case with ground truth clustering
        if(is.null(nc)){
            nc_estim <- length(unique(gt.clust))
            nc <- nc_estim
            message("The number of clusters is automatically set to the number of clusters in the ground truth clustering!")
        } else if(length(nc) == 1){
            if(nc != length(unique(gt.clust))){
                stop("The number of clusters is different from given clustering!")
            } else {
                nc_estim <- nc
                message("The number of clusters is set to ", nc, "!")
            }
        } else if(length(nc) > 1){
            if(length(unique(gt.clust)) %in% nc){
                if(int_method == "SNF"){
                  nc_estim <- SNFtool::estimateNumberOfClustersGivenGraph(W_int, NUMC=nc)$`Eigen-gap best`
                } else if(int_method == "NEMO"){
                  nc_estim <- NEMO::nemo.num.clusters(W_int, NUMC=nc)
                }
                warning("The number of clusters is automatically estimated in range: ", min(nc), ":", max(nc))
            } else {
                stop("The number of clusters to evaluate does not include the number of clusters in the ground truth clustering!")
            }
        }
        
    }
    
    return(list(nc_estim=nc_estim, nc_range=nc))
}


#' Rank features by NMI (modified from SNFtool)
#' 
#' @description
#' Calculates the normalized mutual information (NMI) score between each 
#' features clustering and the clustering of the fused matrix W. Each feature
#' is ranked based on how similar it is to the clustering of the fused matrix.
#'
#' This function is modified from SNFtool::rankFeaturesByNMI to take as input 
#' a precomputed clustering of the fused matrix W and compute the scaled 
#' exponential euclidean kernel as recommended in SNFtool for continuous values.
#'
#' @param data list. A list of matrices.
#' @param clustering vector. Clustering of the fused matrix W.
#' @param K vector. Vector with the number of neighbours for each omic.
#' @param sigma integer. Variance for local model.
#'
#' @returns A list containing NMI score for each feature from all data types
#' and their NMI score ranks.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#'
rankFeaturesByNMI_mod <- function(data, clustering, K, sigma=0.5){  
    
    stopifnot(class(data) == "list")
    
    NUM.OF.DATA.TYPES <- length(data)
    NMI.scores <- vector(mode="list", length=NUM.OF.DATA.TYPES)
    NMI.ranks <- vector(mode="list", length=NUM.OF.DATA.TYPES)
    #num.of.clusters.fused <- estimateNumberOfClustersGivenGraph(W)[[1]]
    num.of.clusters.fused <- length(unique(clustering)) # [JG]
    #clustering.fused <- spectralClustering(W, num.of.clusters.fused)
    clustering.fused <- clustering # [JG]
    
    
    for (data.type.ind in 1:NUM.OF.DATA.TYPES){
        
        NUM.OF.FEATURES <- dim(data[[data.type.ind]])[2] 
        NMI.scores[[data.type.ind]] <- vector(mode="numeric", 
                                              length=NUM.OF.FEATURES)
        
        for (feature.ind in 1:NUM.OF.FEATURES){
            #affinity.matrix <- affinityMatrix(
            #    dist2(as.matrix(data[[data.type.ind]][, feature.ind]), 
            #          as.matrix(data[[data.type.ind]][, feature.ind])))
            
            affinity.matrix <- affinityMatrix(
                (SNFtool::dist2(as.matrix(data[[data.type.ind]][, feature.ind]), 
                      as.matrix(data[[data.type.ind]][, feature.ind])))^(1/2), K=K[data.type.ind], sigma=sigma#[JG]
                )    
            
            clustering.single.feature <- spectralClustering(affinity.matrix, 
                                                            num.of.clusters.fused)
            
            # [JG]: Subset in case the omic has less samples than the fused matrix
            clustering.fused_sub <- clustering.fused[rownames(affinity.matrix)]
            
            NMI.scores[[data.type.ind]][feature.ind] <- calNMI(clustering.fused_sub, 
                                                               clustering.single.feature)      
        }    
        
        NMI.ranks[[data.type.ind]] <- rank(-NMI.scores[[data.type.ind]],
                                           ties.method="first")
    }
    
    return(list(NMI.scores, NMI.ranks))  
}

#' Rank features by NMI
#'
#' @param mat_data list. List of dataframes (samples x features), one for each modality.
#' @param clustering vector. Clustering of the fused matrix W.
#' @param K vector. Vector with the number of neighbours for each omic.
#' @param sigma integer. Variance for local model.
#'
#' @return Dataframe in long-format with score and ranking for each features in each modality.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
compute_feature_importance <- function(mat_data, clustering, K, sigma=0.5){
    
    feature_ranks <- rankFeaturesByNMI_mod(mat_data, clustering, K=K, sigma=sigma)
    
    res <- list()
    for(i in 1:length(mat_data)){
        res[[i]]  <- as.data.frame(cbind(feature_ranks[[1]][[i]], feature_ranks[[2]][[i]]))
        colnames(res[[i]]) <- c("score", "rank")
        res[[i]]$modality <- names(mat_data)[i]
        res[[i]]$features <- colnames(mat_data[[i]])
        
        res[[i]] <- res[[i]] %>% dplyr::arrange(rank)
        
    }
    
    res <- do.call(rbind, res)
    
    return(res)
}


#' Clinical variables enrichment and survival analysis
#' 
#' @description
#' This function performs:
#' - Enrichment analysis for clinical variables using Chi-squared test for 
#' categorical variables and Kruskal-Wallis test for numerical variables
#' - Survival analysis using the log-rank test
#' 
#' The p-values for enrichment of clinical variables are corrected for multiple
#' testing using Benjamini-Hochberg procedure.
#' A Kaplan-Meier (KM) plot us saved for each endpoint and number of clusters considered.
#' 
#'
#' @param clustering dataframe. Dataframe (samples x number of clusters).
#' @param metadata dataframe. Dataframe with metadata (samples x variables). 
#' @param enrich_vars vector. Names of the variables for enrichment analysis.
#' @param surv_vars vector. Names of the endpoint (event and time) for 
#' survival analysis.
#' @param file_path string. Path to save results
#'
#' @return KM plots and p-values of statistical tests.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#'
enrich_surv_analysis <- function(clustering, metadata, enrich_vars=c(), 
                                 surv_vars=c(), file_path){
    # Compute enrichment analysis and log-rank test for each cluster
    # with respect to the selected variables
    
    res <- list(enrich_res=NULL, surv_res=NULL)
    
    ###### Compute enrichment analysis
    if(length(enrich_vars) > 0){
        
        enrich_res <- apply(clustering, 2, function(cl) {
            message("Computing enrichment analysis for cluster K: ", length(unique(cl)))
            ps <- sapply(enrich_vars, function(v) {
                
                if(is.character(metadata[[v]])){
                    print(cl)
                    enrich_data <- data.frame(var=metadata[[v]], clust=as.factor(cl))
                    enrich_data
                    enrich_data <- enrich_data[complete.cases(enrich_data), ] # Consider only complete.cases
                    
                    
                    
                    # apply chi-squared test
                    #res.sim <- chisq.test(as.factor(metadata[,v]), as.factor(cl), 
                    #                    simulate.p.value = T, B=10000)
                    message("\tChi-squared test on variable ", v)
                    p.val <- chisq.test(enrich_data$var, enrich_data$clust, 
                                        simulate.p.value = F)$p.value
                    
                } else if(is.numeric(metadata[[v]])){
                    
                    enrich_data <- data.frame(var=metadata[[v]], clust=as.factor(cl))
                    enrich_data <- enrich_data[complete.cases(enrich_data), ] # Consider only complete.cases
                    
                    # apply kruskal-wallis test
                    message("\tKruskal-Wallis test on variable ", v)
                    p.val <- kruskal.test(enrich_data$var, enrich_data$clust)$p.value
                }
                
            })
            
            # Adjust p-values
            print(ps)
            ps.adj = p.adjust(ps, method = "BH")
        })
        
        res$enrich_res <- enrich_res
    }
    
    ####### Compute survival analysis
    message("\n")
    if(length(surv_vars) > 0){
        #Get names of endpoint to analyze
        endpoint_names <- surv_vars
        
        remove_suffix <- function(x) {
            gsub("\\.time|\\.event", "", x)
        }
        
        endpoint_names <- unique(remove_suffix(endpoint_names))
        
        # Compute log-rank test
        # event is given as input as yes and no
        
        surv_res <- apply(clustering, 2, function(cl) {
            message("Computing survival analysis for cluster K: ", length(unique(cl)))
            
            ps <- sapply(endpoint_names, function(v) {
                
                # plot KM curves
                message("\tComputing log-rank test for endpoint ", v)
                
                surv_data <- metadata[, c(paste0(v, ".time"), paste0(v, ".event"))]
                surv_data$clust <- cl
                surv_data <- surv_data[complete.cases(surv_data),] # Consider only complete.cases
                
                if(any(surv_data[, paste0(v, ".event")] %in% c('yes', 'no'))){
                    surv_data[, paste0(v, ".event")] = as.numeric(surv_data[, paste0(v, ".event")] == 'yes')
                }
                
                formula.text <- paste0('Surv(', paste0(v, ".time,"),paste0(v, ".event"),') ~ clust')
                form <- as.formula(formula.text)
                
                km <- surv_fit(formula=form, data=surv_data)
                
                ggsurvplot(km, data=surv_data, title= v, conf.int = T, ggtheme = theme_gray()) 
                ggsave(paste0("KM_", v, "_K",length(unique(cl)), ".png"),  bg="white", 
                       path=file_path, width=15, height=10, units="cm")
                
                p.val <- survdiff(formula=form, data=surv_data)$pvalue
            })
        })
        
        res$surv_res <- surv_res
    }
    
    return(res)
}


#' Plot internal quality indices
#'
#' @description
#' 
#'  Plot the Silhouette score and Dunn index w.r.t. the number of clusters
#'  evaluated. The plot is saved in the directory "clustering_metrics", which
#'  is automatically created in the current position if not exists.
#'
#' @param int.val.idx dataframe. Dataframe (indices x number of clusters).
#' @param nc_range vector. Vectors with the number of clusters evaluated.
#'
#' @return Save a png plot of the Silhouette score, Dunn index and Entropy.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#'
plot_int.val.idx <- function(int.val.idx, nc_range){
    
    to_plot <- int.val.idx[c("avg.silwidth", "dunn", "entropy", "ch"), , drop=F] %>% 
        tibble::rownames_to_column(var="metrics") %>% 
        tidyr::pivot_longer(!metrics, names_to="nc") %>%
        dplyr::mutate_at('value', as.numeric) %>% 
        dplyr::mutate_at(c("metrics", "nc"), as.factor) %>%
        dplyr::mutate(metrics = dplyr::case_when(
            metrics == "avg.silwidth" ~ "Mean Silhouette Coefficient",
            metrics == "dunn" ~ "Dunn Index",
            metrics == "entropy" ~ "Entropy",
            metrics == "ch" ~ "Calinski-Harabasz Index",
            TRUE ~ metrics  # Keep other values unchanged
        )) %>%
        dplyr::rename("Indices" = "metrics")
        
    
    #to_plot$K <- factor(to_plot$K, levels=paste0("nc", nc_range))
    
    dir.create("clustering_metrics", showWarnings = FALSE)
    
    ggplot(to_plot, aes(x=nc, y=value, color=Indices, group=Indices)) +
        geom_point() +
        geom_line() +
        theme_minimal() +
        facet_wrap(~Indices, scales="free_y") +
        labs(title="Internal quality indices",
             x="Number of clusters", y="") +
        theme(legend.position="bottom", plot.title = element_text(hjust = 0.5))
    ggsave(file.path("clustering_metrics", "internal_quality_indices.png"), 
           width=10, height=5, bg="white")
    
    # Save data
    write.csv(pivot_wider(to_plot, names_from = nc), 
              file="clustering_metrics/internal_quality_indices.csv", row.names=T)
    
}


#' Plot external quality indices
#'
#' @description
#' 
#'  The plot is saved in the directory "clustering_metrics", which
#'  is automatically created in the current position if not exists.
#'
#' @param ext.val.idx dataframe. Dataframe (indices x number of clusters)
#' @param nc_range vector. Vectors with the number of clusters evaluated.
#' @param types vector. Vector with the external quality indices to plot.
#'
#' @return Save a png plot of external quality indices.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
#'
plot_ext.val.idx <- function(ext.val.idx, nc_range, 
                             types = c("ARI", "AMI", "NVI", "NMI", "V-measure")){
    
    ext.val.idx <- ext.val.idx %>% 
        tidyr::pivot_longer(all_of(types), names_to = "Indices", values_to = "scores") %>%
        dplyr::mutate(dplyr::across(where(is.character), as.factor))
    
    ggplot(ext.val.idx, aes(x=nc, y=scores, color=Indices, group=Indices)) +
        geom_point() +
        geom_line() +
        theme_minimal() +
        facet_wrap(~Indices, scales="free_y") +
        labs(title="External quality indices",
             x="Number of clusters", y="") +
        theme(legend.position="bottom", plot.title = element_text(hjust = 0.5))
    
    ggsave(file.path("clustering_metrics", "external_quality_indices.png"), 
           width=10, height=5, bg="white")
    
    # Save data
    write.csv(pivot_wider(ext.val.idx, names_from = nc, values_from = scores), 
              file="clustering_metrics/external_quality_indices.csv", row.names=F)

}


#' Plot alluvial diagram
#'
#' @description Alluvial plot to follow the flow of samples between predicted clustering, 
#' ground truth clustering and variable of interest. In particular, one plot
#' is automatically saved for each variable of interest.
#' 
#' If no ground truth clustering or variable of interest is provided, no plot 
#' is generated and a message is printed.
#' 
#'
#' @param clust_df dataframe. Dataframe containing the cluster assignments.
#' @param gt.clust factor. Factor containing the ground truth clustering (can be NULL).
#' @param var_explain dataframe. Dataframe containing the variables to explain 
#' the clustering (can be NULL but no plot will be created).
#' 
#' @param save_path string. Path to save the plot.
#'
#' @return Save a png plot of the alluvial diagram.
#' 
#' @author Jessica Gliozzo
#' 
#' @export


#PATH FOR GGPLOT2
#fill_alpha() function is no longer available in new ggplot versions
if (!exists("fill_alpha")) {
  fill_alpha <- function(fill, alpha) {
    scales::alpha(fill, alpha)
  }
}

plot.alluvial <- function(clust_df, gt.clust, var_explain, save_path="."){
    
    # Check if gt.clust or var_explain are NULL
    if(is.null(gt.clust) & is.null(var_explain)){
        message("Both gt.clust and var_explain are NULL. No alluvial plot will be generated.")
        return(NULL)
    } else if(is.null(var_explain)){
        message("var_explain is NULL. No alluvial plot will be generated.")
        return(NULL)
    }
    
    # check if var_explain has columns
    if(ncol(var_explain) == 0){
      message("var_explain is empty. No alluvial plot will be generated.")
      return(NULL)
    }
  
    # Check that clust_df, gt.clust and var_explain have the same order
    sample_names <- list(names(gt.clust), rownames(var_explain))
    sample_names <- sample_names[ !sapply(sample_names, is.null) ]
    same_order <-  all(vapply(sample_names, function(x) identical(x, rownames(clust_df)), logical(1)))
    if(!same_order){
      stop("The rownames of clust_df, gt.clust and var_explain must be the same!")
    }
  
    names(clust_df) <- "Pred. Clusters"
    
    # detect axes and var
    axes = colnames(clust_df)
    vars = colnames(var_explain)
    
    # add var to clust_df
    if(!is.null(gt.clust)){
        clust_df = cbind(clust_df, gt.clust, var_explain)
        axes = c(axes, "Clusters")
    } else {
        clust_df = cbind(clust_df, var_explain)
    } 
    
    # Fix columns type
    clust_df = clust_df %>% 
      { if ("gt.clust" %in% names(.)) rename(., "Clusters" = "gt.clust") else . } %>%
         dplyr::mutate(across(where(is.character), as.factor)) %>%
         dplyr::mutate_at(axes, as.factor)
    
    # loop over vars
    for(v in vars){
        
        # Bin continuous variables
        if(is.numeric(clust_df[[v]])){
            clust_df[[v]] = cut(clust_df[[v]], 
                                breaks =  unique(pretty(clust_df[[v]])), 
                                include.lowest = T)
        }
        
        # get frequency of each combination of cluster assignments
        df = clust_df[, c(axes, v)] %>% dplyr::group_by_all() %>% dplyr::count(name='Frequency')
        
        # Create aes string 
        aes_text = paste0("aes(", paste0(paste0("axis", 1:length(axes), " = .data[[axes[", 1:length(axes), "]]]"), collapse = ", "), ", ",
                          paste0("axis", length(axes)+1, " = ", v), ", y = Frequency)")
        
        # plot
        ggplot(data = df,
               # aes(axis1 = .data[[axes[1]]], axis2 = .data[[axes[2]]], axis3 = .data[[axes[3]]], 
               #     y = freq)) +
               eval(parse(text = aes_text))) +
            scale_x_discrete(limits = c(axes, v), expand = c(.2, .05)) +
            #xlab("Clusters and explain variable") +
            geom_alluvium(aes(fill = !!sym(v))) +
            geom_stratum() +
            #geom_text(stat = "stratum", aes(label = after_stat(stratum)), size= 2.5) +
            ggfittext::geom_fit_text(stat = "stratum", aes(label = after_stat(stratum)), 
                                     width = 1/3, min.size = 4) +
            theme_minimal() +
            ggtitle("Alluvial plot of cluster assignments") + 
            theme(plot.title = element_text(hjust = 0.5))
        
        #save plot
        ggsave(file.path(save_path, paste0("alluvial_plot_", v, ".png")), 
               width = 15, height = 10, units = "cm", dpi = 300, bg = "white")
        
    }
    
}


#' Check rownames and colnames of similarity matrices and metadata
#'
#' @param sim_mat list. List of affinity matrices.
#' @param metadata list. List of metadata dataframes.
#' @param int_method string. Method used for data fusion.
#'
#' @returns A message if the rownames and colnames of the data and metadata are 
#' identical, otherwise it stops.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
check_names <- function(sim_mat, metadata, int_method){
    
    for(omic in names(sim_mat)){
        if(identical(rownames(sim_mat[[omic]]), metadata$ID) & 
           identical(colnames(sim_mat[[omic]]), metadata$ID)){
            message("The rownames and colnames of ", omic, " are identical in data and metadata")
        } else {
            stop("The rownames and colnames of ", omic, " are not identical in data and metadata")
        }
        
    }
    
    # Check if the rownames of different omics are identical
    for(i in 1:(length(sim_mat)-1)){
        if(identical(rownames(sim_mat[[i]]), rownames(sim_mat[[i+1]]))){
            message("\nThe rownames of ", names(sim_mat)[i], " and ", names(sim_mat)[i+1], " are identical")
        } else {
            message("\nThe rownames of ", names(sim_mat)[i], " and ", names(sim_mat)[i+1], " are not identical")
            if(int_method=="SNF"){
              stop("SNF requires the same samples order for all omics. Please check your data.")
            } else if(int_method=="NEMO"){
              message("NEMO does not require the same samples order. Do not worry!")
            }
        }
    }
}


#' Plot heatmap of feature importance
#'
#' @description
#' This function plots a heatmap of the top features. 
#' The features are ordered by their importance score and the samples are
#' ordered by their ground truth clustering (if available) or predicted clustering.
#' 
#'
#' @param data list. List of dataframes (samples x features), one for each modality.
#' @param feat_imp dataframe. Dataframe containing the feature importance scores.
#' @param pred.clust vector. Predicted clustering of the samples.
#' @param gt.clust factor. Factor containing the ground truth clustering (can be NULL).
#' @param top_feat integer. Number of top features to plot for each omic.
#' @param save_path string. Path to save the plot.
#'
#' @returns A heatmap of the top features for each omic, saved as a PNG file.
#' 
#' @author Jessica Gliozzo
#' 
#' @export
plot_imp_heatmap <- function(data, feat_imp, pred.clust, gt.clust=NULL, top_feat=10, 
                               save_path="."){
    
    # To set annotation colors
    pred_cols <- setNames(
      brewer.pal(n = length(levels(pred.clust)), name = "Set1"),
      levels(pred.clust)
    )
    pred_cols <- pred_cols[levels(pred.clust)] #remove unwanted colors when 2 clust, min 3 colors
    
    if(!is.null(gt.clust)){
      gt_cols <- setNames(
        brewer.pal(n = length(levels(gt.clust)), name = "Set2"),
        levels(gt.clust)
      )
      gt_cols <- gt_cols[levels(gt.clust)]  #remove unwanted colors when 2 clust, min 3 colors
    }
  
    # To save results
    all_data <- list()
    all_imp <- list()
    top_anno <- list()
    
    # copy of pred.clust and gt.clust
    pred.clust_orig <- pred.clust
    gt.clust_orig <- gt.clust
    
    # Iterate over omics
    for(omic in names(data)){
        
        data_mat <- data[[omic]]
        
        # subset pred.clust and gt.clust to the samples in data_mat
        if(!is.null(gt.clust)){
            pred.clust <- pred.clust_orig[rownames(data_mat)]
            gt.clust <- gt.clust_orig[rownames(data_mat)]
        } else {
            pred.clust <- pred.clust_orig[rownames(data_mat)]
        }
        
        # extract feature importance for omic
        feat_imp_omic <- feat_imp[feat_imp$modality == omic, ]
        
        # select the top features for each omic
        top_feat_omic <- feat_imp_omic[1:top_feat, ]
        
        # select the data for the top features
        data_mat <- data_mat[, match(top_feat_omic$features, colnames(data_mat))]
        
        # transpose omic data (features x samples)
        data_mat <- t(data_mat)
        
        if(!is.null(gt.clust)){
            
            # reorder gt.clust, pred.clust and data_mat based on gt.clust
            gt.clust <- sort(gt.clust)
            ordered_samples <- names(gt.clust)
            pred.clust <- as.factor(pred.clust[ordered_samples])
            data_mat <- data_mat[, ordered_samples]
            
            # Check if rownames of data are in the same order of gt.clust and pred.clust
            if(identical(colnames(data_mat), names(gt.clust)) & identical(names(gt.clust), names(pred.clust))){
                message(paste(omic, "data, predicted clustering and ground truth clustering have same order!"))
            } else {
                stop(paste(omic, "data, predicted clustering and ground truth clustering have NOT the same order!"))
            }
            
            # store results
            all_data[[omic]] <- data_mat
            all_imp[[omic]] <- top_feat_omic
            
            # Annotations
            df <- data.frame("Pred.Clusters"=pred.clust, "Clusters"=gt.clust)
            top_anno[[omic]] <- HeatmapAnnotation(df=df, col = list("Pred.Clusters" = pred_cols, 
                                                                    "Clusters" = gt_cols))
            
        } else{
            
            # reorder pred.clust and data_mat based on pred.clust
            pred.clust <- as.factor(sort(pred.clust))
            ordered_samples <- names(pred.clust)
            data_mat <- data_mat[, ordered_samples]
            
            # Check if rownames of data are in the same order of  pred.clust
            if(identical(colnames(data_mat), names(pred.clust))){
                message(paste(omic, "data and predicted clustering have same order!"))
            } else {
                stop(paste(omic, "data and predicted clustering have NOT the same order!"))
            }
            
            # store results
            all_data[[omic]] <- data_mat
            all_imp[[omic]] <- top_feat_omic
            
            # Annotations
            df <- data.frame("Pred.Clusters"=pred.clust)
            top_anno[[omic]] <- HeatmapAnnotation(df=df, col = list("Pred.Clusters" = pred_cols))
            
        }
    } # end iteration across omics
    
    
    # make plot
    shared_col_fun <- colorRamp2(c(0, 1), c("white", "green"))
    
    heatmap_list <- lapply(seq_along(all_data), function(i){
        
        right_anno <- rowAnnotation("Score"=all_imp[[i]]$score, 
                                    col = list("Score"=shared_col_fun), 
                                    show_annotation_name = FALSE)
        
        ComplexHeatmap::Heatmap(all_data[[i]], cluster_rows=F, cluster_columns = F,
                                show_column_names = F, name="Values", top_annotation = top_anno[[i]],
                                row_names_side = "left", row_names_gp = gpar(fontsize = 8),
                                right_annotation = right_anno)
    })

    # plot_heatmaps
    for(i in 1:length(heatmap_list)) {
        png(file.path(save_path, paste0(names(data)[i], "_top_features_heatmap.png")), width = 1600, height = 2000, res = 150)
        draw(heatmap_list[[i]],  merge_legend = T, annotation_legend_side = "right")
        dev.off()
    }
    
}



#' Glue function handed to clusterboot for spectral clustering
#'
#' @param dist_mat matrix. Distance matrix (samples x samples), already normalized.
#' @param k integer. Number of clusters.
#'
#' @returns List containing:
#' - result: clustering result
#' - partition: vector of cluster assignments
#' - clusterlist: list of logical membership vectors
#' - nc: number of clusters
#' - nccl: number of proper clusters
#' - clustermethod: clustering method used
#' 
#' @author Jessica Gliozzo
#' 
#' @export
fpc_spectralClustering <- function(dist_mat, k) {
  
  W <- 1 - as.matrix(dist_mat)
  
  part <- SNFtool::spectralClustering(W, K = k)   # SNFtool
  
  # turn vector into list of logical membership vectors, as clusterboot expects
  clist <- lapply(seq_len(k), function(j) part == j)
  
  list(result = part,
       partition = part,
       clusterlist = clist,
       nc = k,           # number of clusters incl. “noise” (none here)
       nccl = k,         # number of “proper” clusters
       clustermethod = "spectralClustering")
}


#' Extract subset statistics from clusterboot results
#'
#' @param cb_list list. List of clusterboot results.
#' @param k_vec vector. Vector with the number of clusters evaluated.
#'
#' @returns Dataframe with the subset statistics for each number of clusters:
#' - subsetmean: clusterwise means of Jaccard similarities
#' - subsetsd: clusterwise standard deviations of Jaccard similarities
#' - subsetbrd: clusterwise number of times a cluster has been dissolved
#' - subsetrecover: clusterwise number of times a cluster has been recovered
#' 
#' @author Jessica Gliozzo
#' 
#' @export
extract_subset_stats <- function(cb_list, k_vec) {
  
  out <- lapply(seq_along(cb_list), function(i) {
    cb <- cb_list[[i]]
    
    data.frame(
      k              = k_vec[i],
      subsetmean     = cb$subsetmean, 
      subsetsd       = apply(cb$subsetresult, 1, sd),
      subsetbrd      = cb$subsetbrd,  
      subsetrecover  = cb$subsetrecover
    )
  })
  
  res <- do.call(rbind, out)
  
  return(res)
}


#' @title NEMO affinity graph (fix bug)
#' @name nemo.affinity.graph_fix
#' @description Constructs a single affinity graph measuring similarity across different omics.
#' This function is modified from NEMO::nemo.affinity.graph to fix a bug in the original implementation 
#' that does not allow to provide as k a vector of numbers for each omic.
#'
#' See paper for details:
#' Nimrod Rappoport, Ron Shamir, NEMO: cancer subtyping by integration of partial 
#' multi-omic data, Bioinformatics, Volume 35, Issue 18, September 2019, Pages 
#' 3348–3356, https://doi.org/10.1093/bioinformatics/btz058
#' 
#' @param raw.data A list of the data to be clustered, where each an entry is a dataframe of features x samples.
#' @param k The number of neighbors to use for each omic. It can either be a number, a list of numbers
#' or NA. If it is a number, this is the number of neighbors used for all omics. If this is a list,
#' the number of neighbors are taken for each omic from that list. If it is NA, each omic chooses the
#' number of neighbors to be the number of samples divided by NUM.NEIGHBORS.RATIO.
#' @return A single matrix measuring similarity between the samples across all omics.
#' 
#' @author Nimrod Rappoport, Ron Shamir
#' 
#' @export
nemo.affinity.graph_fix <- function(raw.data, k=NA) {
  if (any(is.na(k))) {
    k = as.numeric(lapply(1:length(raw.data), function(i) round(ncol(raw.data[[i]]) / NUM.NEIGHBORS.RATIO)))
  } else if (length(k) == 1) {
    k = rep(k, length(raw.data))
  }
  sim.data = lapply(1:length(raw.data), function(i) {affinityMatrix(dist2(as.matrix(raw.data[[i]]),
                                                                          as.matrix(raw.data[[i]])), k[i], 0.5)})
  affinity.per.omic = lapply(1:length(raw.data), function(i) {
    sim.datum = sim.data[[i]]
    non.sym.knn = apply(sim.datum, 1, function(sim.row) {
      returned.row = sim.row
      threshold = sort(sim.row, decreasing = T)[k[i]]
      returned.row[sim.row < threshold] = 0
      row.sum = sum(returned.row)
      returned.row[sim.row >= threshold] = returned.row[sim.row >= threshold] / row.sum
      return(returned.row)
    })
    sym.knn = non.sym.knn + t(non.sym.knn)
    return(sym.knn)
  })
  patient.names = Reduce(union, lapply(raw.data, rownames))
  num.patients = length(patient.names)
  returned.affinity.matrix = matrix(0, ncol = num.patients, nrow=num.patients)
  rownames(returned.affinity.matrix) = patient.names
  colnames(returned.affinity.matrix) = patient.names
  
  shared.omic.count = matrix(0, ncol = num.patients, nrow=num.patients)
  rownames(shared.omic.count) = patient.names
  colnames(shared.omic.count) = patient.names
  
  for (j in 1:length(raw.data)) {
    curr.omic.patients = rownames(raw.data[[j]])
    returned.affinity.matrix[curr.omic.patients, curr.omic.patients] = returned.affinity.matrix[curr.omic.patients, curr.omic.patients] + affinity.per.omic[[j]][curr.omic.patients, curr.omic.patients]
    shared.omic.count[curr.omic.patients, curr.omic.patients] = shared.omic.count[curr.omic.patients, curr.omic.patients] + 1
  }
  
  final.ret = returned.affinity.matrix / shared.omic.count
  lower.tri.ret = final.ret[lower.tri(final.ret)]
  final.ret[shared.omic.count == 0] = mean(lower.tri.ret[!is.na(lower.tri.ret)])
  
  return(final.ret)
}
