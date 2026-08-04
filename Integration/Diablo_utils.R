# Functions used in main script "SNF_int.R"
# Author: Cristian Iperi


#' The input for the integration analyses generated in the single omics pipeline
#' are transformed in this function of obtain a consensus format. 
#'
#' @param matrix matrix format. It contains the matrix of data including the first column ID.
#' The ID column represents the variables in transcriptomics and methylomics data, while the sample ID
#' in metabolomics and undefined omics (the differences in the format are due to upstream constraint)
#'
#' @return A list with including lists of two element (one per omics input):
#' - [1]: a the data in a dataframe format for each omics
#' - [2]: Variable Identifier in a vector format for each omics.
#' 
#' @author Cristian Iperi
#' 
#' @export
#### UNDEFINED FUNCTION

Undefined_processing <-function(matrix,mer){
  x<- as.data.frame(colnames(matrix))
  
  matrix <- matrix %>% filter(CONDITION == args[2] | CONDITION == args[1]) #Args[1] = Control condition name/  args[2] = Disease condition name
  
  list <- list()        
  
  sample<-matrix$ID
  matrix<-matrix %>% dplyr::select(!c(CONDITION, ID))
  #colnames(matrix) <- gsub("peak","peakS",colnames(matrix))
  list[[1]] <- as.data.frame(matrix)
  rownames(list[[1]]) <- sample
  
  sample_undefined_EASY <- colnames(matrix)
  list[[2]] <- sample_undefined_EASY
  
  
  return(list)
  
}




#### METABOLOMICS FUNCTION

Metabolomics_processing <-function(annotation,matrix,mer){
  #matrix[,-1:-2]<- apply(matrix[,-1:-2],2,log)
  annotation <- annotation %>% distinct(NAME.x, .keep_all = TRUE)
  annotation$NAMES <- paste(annotation$NAME.x, annotation$Name,sep = "/")
  x<- as.data.frame(colnames(matrix))
  x$id  <- 1:nrow(x)
  colnames(x)[1] <- "Peaks"
  xx<-merge(x,annotation,by.x="Peaks", by.y ="NAME.x", all.x = TRUE)
  xx <- xx[order(xx$id), ]
  
  
  for (i in 1:nrow(xx)){
    if (!is.na(xx$NAMES[i])){
      xx$Peaks[i]<- xx$NAMES[i]
    }
  }
  
  serum_metabolomics_EASY <- xx$Peaks
  
  list <- list(annotation)
  
  
  matrix <- matrix %>% filter(CONDITION == args[2] | CONDITION == args[1]) 
  sample<-matrix$ID
  matrix<-matrix %>% dplyr::select(!c(CONDITION, ID))
  #colnames(matrix) <- gsub("peak","peakS",colnames(matrix))
  list[[1]] <- as.data.frame(matrix)
  rownames(list[[1]]) <- sample
  
  #ADD IF STATEMENT IF 2 METABOLOMICS
  #serum_metabolomics$feature <- gsub("peak",paste("peak_",mer, sep = ""), serum_metabolomics$feature)
  
  list[[2]] <- serum_metabolomics_EASY
  
  return(list)
  
}




#### TRANSCRIPTOMIC FUNCTION


Transcriptomics_processing <-function(annotation,matrix,mer){
  list <-list()
  matrix<-merge(matrix, MART, by.x="ID", by.y="Gene.stable.ID")
  matrix<- matrix %>% arrange(desc(variance))
  matrix <-matrix[1:Max_features,]
  
  Bcell_RNAseq_EASY <- paste(matrix$`Gene name`, matrix$`Gene.name`, sep = "/")
  Bcell_RNAseq_VIEW <- as.character(matrix$`Gene name`)
  
  matrix <- t(matrix)
  colnames(matrix) <- matrix[1,]
  matrix <- as.data.frame(matrix)
  matrix <- matrix[c(-1, -(nrow(matrix) - 1),-nrow(matrix)),]
  
  Bcell_RNAseq_VIEW <-make.unique(Bcell_RNAseq_VIEW, sep = "_")
  colnames(matrix) <- Bcell_RNAseq_VIEW
  list[[1]] <- matrix
  list[[2]] <- Bcell_RNAseq_VIEW
  
  return(list)
}



#### METHYLOMICS FUNCTION


Methylomics_processing <-function(annotation,matrix,metadata,mer){
  annotation <- annotation[,c("gene","CpG_island")]
  matrix<-merge(matrix, annotation, by.x="ID", by.y="CpG_island")
  
  matrix <-matrix[1:Max_features,]
  
  Methylome_WB_VIEW <- as.character(matrix$gene)
  Methylome_WB_EASY <- paste(matrix$ID, matrix$gene, sep = "/")
  
  list <-list()
  
  matrix <- t(matrix)
  colnames(matrix) <- matrix["ID",]
  matrix <- as.data.frame(matrix)
  y <-rownames(matrix) %in% c("ID","variance")
  matrix <- matrix[!y,]
  matrix <- matrix[c(-(nrow(matrix) - 1),-nrow(matrix)),]
  
  Methylome_WB_EASY <-make.unique(Methylome_WB_EASY, sep = "_")
  colnames(matrix) <- Methylome_WB_EASY
  list[[1]] <- matrix
  list[[2]] <- Methylome_WB_EASY
  
  return(list)
}

