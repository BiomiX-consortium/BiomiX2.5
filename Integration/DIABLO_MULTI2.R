#use_python("/usr/bin/python3", required=TRUE)

#
# # #MANUAL INPUT
# args = as.list(c("BLymphocytes","SLE"))
# args[1] <-"mutated"
# args[2] <-"unmutated"
# args[3] <-"C:/Users/crist/Desktop/BiomiX2.5"

# directory <-args[3]
# renv::load(paste(directory,"_INSTALL",sep="/"))

#upload libraries

library(data.table)
library(vroom)
library(dplyr)
library(tidyverse)
library(DESeq2)
library(caret)
library(rlist)
library(mixOmics)
library(reticulate)
library(readxl)



MART <- vroom(paste(directory,"/Integration/x_BiomiX_DATABASE/mart_export_37.txt",sep=""), delim = ",")
myList <- list()

COMMAND <- vroom(paste(directory,"COMMANDS.tsv",sep="/"), delim = "\t")
COMMAND_MOFA <- vroom(paste(directory,"COMMANDS_MOFA.tsv",sep="/"), delim = "\t")
COMMAND_ADVANCED <- vroom(paste(directory,"COMMANDS_ADVANCED.tsv",sep="/"), delim = "\t")
Max_features <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[3])
DIR_METADATA <- readLines(paste(directory,"directory.txt",sep="/"))
Cell_type <- as.character(COMMAND_MOFA[1,2])
Contribution_threshold <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_MOFA[3])
n_factor<-as.numeric(COMMAND_MOFA[3,2])
n_iteration <- as.numeric(COMMAND_ADVANCED$ADVANCED_OPTION_DIABLO_OPTIONS[2])
selection_features_DIABLO <- as.character(COMMAND_ADVANCED$ADVANCED_OPTION_DIABLO_OPTIONS[1])
design_input <- as.character(COMMAND_ADVANCED$ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN[1])





if (grepl("\\.xlsx$|\\.xls$", DIR_METADATA)) {
        METADATA <- read_excel(DIR_METADATA)
        print("Metadata Excel File read successfully!")
}else{
        METADATA <-vroom(DIR_METADATA , delim = "\t", col_names = TRUE)}



#### UNDEFINED FUNCTION

Undefined_processing <-function(matrix,mer){
        x<- as.data.frame(colnames(matrix))
        
        matrix <- matrix %>% filter(CONDITION == args[2] | CONDITION == args[1])
        
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
        
        Bcell_RNAseq_EASY <- paste(matrix$`Gene name`, matrix$ID, sep = "/")
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




##### REARRANGEMENT INPUT1 DATA ----

myList <- list()
names_X <- c()
for (i in 1:length(COMMAND$INTEGRATION)){


#i <- 3 #N_input
type <- COMMAND$DATA_TYPE[i]
#Prova ad aggiungere questa linea di codice e a sostituire COMMAND con i possibili
#input sulla base del database di comandi di riferimento. 

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
names_X<- append(COMMAND$LABEL[i],names_X)
}

if(COMMAND$DATA_TYPE[i] == "Transcriptomics"){
        
        print(args[1])
        directory2 <- paste(directory,"/Integration/INPUT/", COMMAND$LABEL[i],"_",args[1],"_vs_", args[2], sep ="")
        Wholeblood_RNAseq <-  vroom(paste(directory2, "/", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], "_normalized_vst_variance.tsv",sep = ""), delim = "\t") #read normalization only
        Wholeblood_metadata <-  vroom(paste(directory2, "/","/Metadata_",COMMAND$LABEL[i], "_", args[1],".tsv",sep = ""), delim = "\t")
        INPUTX<-Transcriptomics_processing(Wholeblood_metadata,Wholeblood_RNAseq,COMMAND$LABEL[i])
        assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
        myList <- list.append(myList,INPUTX[[1]])
        names(myList)
        names_X<- append(COMMAND$LABEL[i],names_X)
        
        
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
        names_X<- append(COMMAND$LABEL[i],names_X)
        
        
}
        
if(COMMAND$DATA_TYPE[i] == "Undefined"){        
        directory2 <- paste(directory,"/Integration/INPUT/", "Undefined_", COMMAND$LABEL[i], "_",args[1],"_vs_", args[2], sep ="")
        samples_undefined <- vroom(paste(directory2,"/Undefined_",COMMAND$LABEL[i], "_MOFA.tsv", sep = ""), delim = "\t")
        INPUTX<-Undefined_processing(samples_undefined,COMMAND$LABEL[i])
        assign(paste("INPUT", i, "_visual", sep=""),INPUTX[[2]])
        myList <- list.append(myList,INPUTX[[1]])
        names_X<- append(COMMAND$LABEL[i],names_X)
        
}        
}
}

#### CHOISE MERGING ----

# Step 1: Convert matrices to data frames
X <- lapply(myList, as.data.frame)

# Step 2: Find common rows
common_rows <- Reduce(function(x, y) intersect(rownames(x), rownames(y)), myList)

X <- lapply(X, function(df) {
        x<-df[rownames(df) %in% common_rows,]
   
})

for (i in 1:length(X)){
        Y<-X[[i]]
        row1 <- rownames(X[[i]])
        col1 <- colnames(X[[i]])
        X[[i]]<- apply(as.matrix(X[[i]]), 2, as.numeric)
        colnames(X[[i]]) <- col1
        rownames(X[[i]]) <- row1
}


METADATA<-METADATA[METADATA$ID %in% common_rows,]
METADATA <- METADATA[order(factor(METADATA$ID, levels = common_rows)), ]
names_X <- rev(names_X)
names(X) <- names_X
Y <- METADATA$CONDITION # use the subtype as the outcome variable

if (design_input != "X"){
  
  if (grepl("\\.xlsx$|\\.xls$", design_input)) {
    design_full <- read_excel(design_input)
    print("Metadata Excel File read successfully!")
  }else{
    design_full <-vroom(design_input , delim = "\t", col_names = TRUE)}
  
  }else{
    
    design_full <- matrix(1, nrow = length(X), ncol = length(X))
    diag(design_full) <- 0
    colnames(design_full) <- rownames(design_full) <- names(X)
  
}



result.diablo <- block.plsda(X, Y, ncomp = n_factor, max.iter = n_iteration, design = design_full ) # run the method



#DOWNSTREAM ANALYSIS

dir.create(path = paste(directory,"/Integration/OUTPUT/","DIABLO", "_", args[1] ,"_vs_", args[2], "_",as.numeric(COMMAND_MOFA[3,2]),"_factors" ,sep="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")

directory2 <- paste(directory,"/Integration/OUTPUT/","DIABLO", "_", args[1] ,"_vs_", args[2], "_",as.numeric(COMMAND_MOFA[3,2]),"_factors" ,sep="") 
setwd(directory2)

pdf(file= paste("DIABLO", "_", args[1] ,"_vs_", args[2],".pdf", sep =""), width = 20, height = 9)


plotIndiv(result.diablo) # plot the samples



#Try to compute the cutoff safely
cutoff_value <- tryCatch({
        p <-plotVar(result.diablo , cutoff = 0.7) # plot the variables
        print(p)

}, error = function(e) {
        # If an error occurs, set a default cutoff (e.g., 0)
        message("Error in cutoff calculation: ", e$message)
        return(0)  # Default cutoff if error happens
})



#Grid plot variance explained


library(ggplot2)
library(tidyr)

total <- as.data.frame(result.diablo[["prop_expl_var"]])
total<-rbind(total,apply(total,2,sum))
rownames(total)[nrow(total)] <- "total"

# Convert 'total' to a tidy data frame
variance_data <- total[,-ncol(total)] %>%
        as.data.frame() %>%
        rownames_to_column(var = "Component") %>%
        pivot_longer(cols = -Component, names_to = "Omics", values_to = "VarianceExplained") %>%
        filter(Component != "total") # Exclude 'total' row if unnecessary

# Calculate limits for color scaling
min_r2 <- min(variance_data$VarianceExplained) * 100
max_r2 <- max(variance_data$VarianceExplained) * 100

# Create a grid plot (heatmap)
p <-ggplot(variance_data, aes(x = Omics, y = Component)) + 
        geom_tile(aes(fill = VarianceExplained * 100), color = "black") +  # Variance as percentage
        scale_fill_gradientn(colors = c("gray97", "darkblue"), limits = c(min_r2, max_r2), name = "R2 (%)") +
        labs(x = "Omics Block", y = "Component", title = "Variance Explained by Omics and Component") +
        theme_classic() +
        theme(
                axis.text = element_text(size = 12),
                axis.line = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_text(size = 14),
                legend.position = "right",
                legend.title = element_text(size = 12),
                legend.text = element_text(size = 10)
        )


print(p)

#SAVING THE SIGNIFICANCY

# Load DIABLO output (assume `result` is your DIABLO model object)
# Extract latent variables (components) for all blocks
variates <- result.diablo$variates

block_info = NULL
for (block in names(variates)){
        block_info<-append(rep(block, ncol(variates$Y)),block_info)
}
block_info <- rev(block_info)

# Combine variates across all blocks (if you want to analyze combined effects)
# Or analyze individual blocks separately
all_components <- do.call(cbind, variates)

# Metadata (replace with your actual metadata)
meta <- METADATA  # Replace `your_metadata` with your metadata object
condition <- METADATA$CONDITION  # Replace `condition` with your grouping variable (e.g., "SLE" and "CTRL")

# Define groups
group1 <- all_components[condition == args[1], ]  # Samples in group 1
group2 <- all_components[condition == args[2], ]  # Samples in group 2

# Initialize results
pval <- NULL
means <- NULL
stdv <- NULL

# Loop over components
for (u in 1:ncol(all_components)) {
        # Wilcoxon test between groups for each component
        res <- wilcox.test(group1[, u], group2[, u], alternative = "two.sided")
        
        # Calculate differences in means and standard deviations
        means <- append(means, mean(group1[, u]) - mean(group2[, u]))
        stdv <- append(stdv, sd(group1[, u]) - sd(group2[, u]))
        
        # Store p-value
        pval <- append(pval, res[["p.value"]])
}

# Adjust p-values for multiple testing
p.adj <- p.adjust(pval, method = "fdr", n = length(pval))

# Create results data frame
results <- data.frame(
        Factors = colnames(all_components),
        p.adj = p.adj,
        means = means,
        standard_deviation = stdv,
        block = block_info
)

# Save results to a file
output_file <- paste("DIABLO_SEPARATION", "_", args[1], "_vs_", args[2], "_factors.tsv", sep = "")
write.table(results, file = output_file, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")





#VIOLIN PLOT

# Load ggplot2 for plotting
library(ggplot2)

# Extract scores (variates) from DIABLO result
variates <- result.diablo$variates

# Create a long-format data frame with scores for each block and component
plot_data <- do.call(rbind, lapply(names(variates), function(block) {
        data.frame(
                Sample = rownames(variates[[block]]),
                Score = as.vector(variates[[block]]),
                Component = rep(colnames(variates[[block]]), each = nrow(variates[[block]])),
                Block = block
        )
}))

# Metadata (replace `METADATA` with your actual metadata object)
plot_data <- merge(plot_data, METADATA, by.x = "Sample", by.y = "ID")

# Create the violin plot
violin_plot <- ggplot(plot_data, aes(x = Component, y = Score, fill = Block)) +
        geom_violin(alpha = 0.7, scale = "width") +
        geom_jitter(aes(color = CONDITION), size = 0.5, width = 0.2, alpha = 0.8) +
        facet_wrap(~ Block, scales = "free_y", nrow = 1) +
        theme_classic() +
        labs(
                title = "Scores per Component and Omics",
                x = "Component",
                y = "Score"
        ) +
        theme(
                strip.text = element_text(size = 12),
                axis.text = element_text(size = 10),
                legend.position = "bottom"
        ) +
        scale_fill_brewer(palette = "Set2") +
        scale_color_brewer(palette = "Set2")

# Display the plot
print(violin_plot)




#TOP CONTRIBUTORS PLOT


# Assuming you already have the result from block.plsda
library(ggplot2)
library(dplyr)
library(tidyr)

# Example: plsda_result <- block.plsda(...)  # Replace with your actual DIABLO result

# Extracting the loadings for each block from the plsda_result
# Let's say plsda_result$loadings contains a list with the loadings for each block

for (val in (1:ncol(result.diablo$loadings$Y))){
loadings_data <- lapply(result.diablo$loadings, function(x) data.frame(Variable = rownames(x), Loading = x[,val]))

# Assign names to the blocks
names(loadings_data) <- names(result.diablo$loadings)


results2<-results[!results$block == "Y",]
results2 <- results2 %>%
        mutate(Factors = gsub("comp", "factor_", Factors))

nice<-which(results2$p.adj < 0.05)
nice<-unique(results2$Factors[nice])
if (length(nice) != 0){
        
        
        pdf(file= paste("Distribution_factors_contributions_", args[1] ,"_vs_", args[2],".pdf", sep =""), width = 10, height = 7)
        
        for (ir in nice) {
                
                y = ir #Factor of interest
                
                print(paste("Significant factor identified:",y,sep=" "))

                #TO IMPLEMENT THE SAMPLE WEIGHT IN EACH FACTOR.
                
                library(tidyr)
                library(dplyr)
                
                # Step 1: Replace "Comp" with "Factor" in the Component column
                plot_data <- plot_data %>%
                        mutate(Component = gsub("comp", "Factor", Component))
                
                # Step 2: Create a new column combining Component and Omics
                df_wide <- as_tibble(plot_data) %>% 
                        filter(Block != "Y") %>% 
                        dplyr::select(Sample, Score, Component, Block) %>% 
                        unite("Factor", Component, Block, sep = "_")
                
                # Step 3: Reshape to wide format
                df_wide <- df_wide %>%
                        pivot_wider(names_from = Factor, values_from = Score) 
                
                # Step 4: Reorder columns so that Factor1_* comes before Factor2_*
                column_order <- c("Sample", sort(setdiff(names(df_wide), "Sample")))
                
                df_wide <- df_wide %>%
                        dplyr::select(all_of(column_order)) %>%
                        dplyr::rename(ID = Sample)
                        
                
                # Print result
                write.table(df_wide, "DIABLO_FACTORS_WEIGHTS.tsv",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
                
                
 
                #Iteration to save the top contributors by 95% quantile and user threshold.
                for (dimen in unique(results2$block)){
                        

                        weights_B <- as.data.frame(loadings_data[[dimen]])
                        a1 <- weights_B
                        colnames(a1)[2] <- "weights"
                        colnames(a1)[1] <- "features"
                        a1$weights<- a1$weights/max(abs(a1$weights))

                        a1$abs_weight<-abs(a1$weights)
                        nov<-quantile(a1$abs_weight, probs=0.95)
                        
                        p <-ggplot(a1, aes(x=abs_weight)) +
                                geom_density(color="darkblue", fill="lightblue")
                        #geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.7, bins=1000) + stat_bin(bins = 1000)
                        
                        p<- p + labs(title=paste( y ," weights in_", dimen, sep="" ))
                        
                        p <-p + theme(text=element_text(size=12,face = "bold"), plot.title = element_text(hjust = 0.5, size=25,face = "bold"))
                        
                        p <- p + geom_vline(aes(xintercept=nov))
                        #p <- p + geom_vline(aes(xintercept=cin))
                        print(p)
                        
                        pos<-a1[a1$abs_weight > nov,]
                        pos
                        pos<-pos[order(pos$abs_weight),]
                        
                        #setwd("/home/cristia/Scrivania/PhD/Bioinfo/MOFA_integration/Databases_copia/Results_Optimization/MOFA_14_factors/Weights_5_percent") 
                        poss<-pos$weights > Contribution_threshold
                        negg<-pos$weights < -Contribution_threshold
                        possis <- as.data.frame(pos[poss,])
                        neggis <- as.data.frame(pos[negg,])
                        possi <- as.data.frame(pos[poss,2])
                        neggi <- as.data.frame(pos[negg,2])
                        
                        write.table(possis, paste(dimen, "_",COMMAND$DATA_TYPE[which(COMMAND$LABEL == dimen)], "_pos_",y,".tsv",sep=""),quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
                        write.table(neggis, paste(dimen, "_",COMMAND$DATA_TYPE[which(COMMAND$LABEL == dimen)], "_neg_",y,".tsv",sep=""),quote = FALSE, row.names = FALSE, col.names = TRUE,sep = "\t")
                        
                        
                }
                
                
        }
        dev.off()
}


# Combine all loadings into a single data frame
combined_loadings <- bind_rows(lapply(names(loadings_data), function(block) {
        data <- loadings_data[[block]]
        data$Omics <- block
        return(data)
}))

# Normalize 'Loading' within each 'Omics' group
combined_loadings <- combined_loadings %>%
        group_by(Omics) %>%
        mutate(weights = Loading / max(abs(Loading), na.rm = TRUE)) %>%
        ungroup()

# Select top N loadings for each block (e.g., top 3 variables)
top_loadings <- combined_loadings %>%
        group_by(Omics) %>%
        top_n(25, abs(Loading)) %>%
        arrange(Omics, desc(abs(Loading)))

# Plotting top loadings using ggplot2
p <-ggplot(top_loadings, aes(x = reorder(Variable, weights), y = weights, fill = Omics)) +
        geom_bar(stat = "identity", show.legend = FALSE) +
        facet_wrap(~ Omics, scales = "free_y") + 
        coord_flip() + 
        labs(title = paste("Top Loadings Factor ", val, " for Each Omics in DIABLO Model"),
             x = "Variable", y = "Loading") +
        theme_minimal() +
        theme(axis.text = element_text(size = 12), 
              axis.title = element_text(size = 14), 
              plot.title = element_text(size = 16, hjust = 0.5))

print(p)
}



setwd(directory2)

line = as.data.frame(result.diablo[["prop_expl_var"]])
line<-round(line,2)
total<-apply(line,2,sum)
total<-rbind(line,total)
rownames(total)[nrow(total)] <- "total"
write.table(total,file=paste("DIABLO_MATRIX", "_", args[1] ,"_vs_", args[2],"_",n_factor,"_factors.tsv", sep =""))




#BAR plot variance explained


library(ggplot2)
library(tidyr)

# Convert to a tidy data frame
total <- total[,-ncol(total)]

total_df <- as.data.frame(total) %>%
        rownames_to_column(var = "Component") %>%
        pivot_longer(cols = -Component, names_to = "Omics", values_to = "VarianceExplained")

# Remove the 'total' row if not needed for the plot
total_df <- total_df[total_df$Component != "total", ]

# Plot the data
p <- ggplot(total_df, aes(x = Component, y = VarianceExplained * 100, fill = Omics)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Variance Explained by Each Factor Across Omics Blocks",
             x = "Component",
             y = "Percentage of Variance Explained (%)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 16, hjust = 0.5),
              axis.title = element_text(size = 14),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 12))



print(p)

dev.off()

