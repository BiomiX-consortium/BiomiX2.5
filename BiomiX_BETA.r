
cat("\n\n\n\n\n          /////////      ///    /////////    ///       ///     ///    ///      ///\n         ///    ///     ///    ///   ///    /////  //////     ///      ///  ///  \n        /////////      ///    ///   ///    ///  ///  ///     ///        ///    \n       ///    ///     ///    ///   ///    ///  ///  ///     ///      ///  ///  \n      ///////////    ///    /////////    ///  ///  ///     ///    ///      ///\n                                                               ///            ///\n                                                                                 /////\n\n\n\n\n\n")

#directory <- "/media/henry/My_Passport/Article_MULTI_BETA"

# #MANUAL INPUT
# args = as.list(c("BLymphocytes","SLE"))
# args[1] <-"C1"
# args[2] <-"C2"
# args[3] <-"C:/Users/crist/Desktop/BiomiX2.5"

old <- Sys.time() # get start time


#Sys.sleep(5)
print("Welcome to BiomiX toolkit")

#Load information about docker
is_in_docker <- as.logical(Sys.getenv("BIOMIX_RUNNING_IN_DOCKER", unset = FALSE))
print(is_in_docker)
#Sys.sleep(5)

print("Information acquired, running the analysis")

#Sys.sleep(5)



#Master function to run the docker systems
generateRunnerFunction <- function(dockerImage, expose_qc_port = FALSE, gpu = FALSE) {  # <-- NUOVO: aggiunto gpu = FALSE
  
  if (!is_in_docker) {
    # --- Caso locale: nessun Docker, comportamento identico a prima ---
    runner <- function(
    command, args = character(), stdout = "", stderr = "", stdin = "",
    input = NULL, env = character(), wait = TRUE,
    minimized = FALSE, invisible = TRUE, timeout = 0
    ) {
      system2(
        command = command, args = args,
        stdout = stdout, stderr = stderr, stdin = stdin, input = input,
        env = env, wait = wait, timeout = timeout
      )
    }
    return(runner)
  }
  
  # --- Caso Docker: delega al container fratello corrispondente ---
  host_shared_path <- Sys.getenv("BIOMIX_HOST_SHARED_PATH")
  if (nchar(host_shared_path) == 0) {
    stop("BIOMIX_HOST_SHARED_PATH non impostata: impossibile montare la ",
         "cartella condivisa nel container di analisi.")
  }
  
  runner <- function(
    command, args = character(), stdout = "", stderr = "", stdin = "",
    input = NULL, env = character(), wait = TRUE,
    minimized = FALSE, invisible = TRUE, timeout = 0
  ) {
    
    # Template: docker run --rm -v <host_path>:/shared [--p 3839:3839] [--gpus all] --entrypoint Rscript
    #           <dockerImage> COMMAND ARGS
    
    port_args <- if (expose_qc_port) c("-p", "3839:3839") else character()
    gpu_args  <- if (gpu) c("--gpus", "all") else character()  # <-- NUOVO
    
    system2(
      command = "docker",
      args = c(
        "run", "--rm",
        "-v", paste0(host_shared_path, ":/shared"),
        port_args,                        # esposta solo se richiesto (QC preview)
        gpu_args,                         # <-- NUOVO — esposta solo se richiesto (accesso GPU)
        "--entrypoint", command,          # es. "Rscript" — sovrascrive l'ENTRYPOINT dell'immagine
        dockerImage,
        shQuote(args)
      ),
      stdout = stdout, stderr = stderr, stdin = stdin, input = input,
      env = env, wait = wait, timeout = timeout
    )
    
  }
  
  return(runner)
}




# --- Rileva se una GPU NVIDIA è realmente passabile al container -----------
gpu_available <- function() {
  res <- tryCatch(
    system2("docker", c("run", "--rm", "--gpus", "all",
                        "nvidia/cuda:12.4.1-base-ubuntu22.04", "nvidia-smi"),
            stdout = FALSE, stderr = FALSE),
    error = function(e) 1L
  )
  identical(res, 0L)
}


# --- Prepara una matrice in formato "PRAMIGO-ready" a partire dal formato
# usato da MOFA/DIABLO, per le omiche che sono trasposte + hanno la colonna
# condizione da eliminare (Metabolomics, Undefined) --------------------------

# --- Prepare a "PRAMIGO-ready" matrix from the MOFA/DIABLO format, for
# omics that are transposed and carry a condition column to drop
# (Metabolomics, Undefined) -------------------------------------------------

prepare_pramigo_matrix <- function(original_path, condition_col_index = 2) {
  
  cat("\n[PRAMIGO prep] Reading original matrix: ", original_path, "\n", sep = "")
  
  df <- read.delim(original_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  cat("[PRAMIGO prep] Original shape: ", nrow(df), " rows x ", ncol(df),
      " columns (samples x [ID + condition + features])\n", sep = "")
  cat("[PRAMIGO prep] First 3 column names: ",
      paste(head(colnames(df), 3), collapse = ", "), "\n", sep = "")
  cat("[PRAMIGO prep] Dropping column ", condition_col_index,
      " (condition) - name was: '", colnames(df)[condition_col_index], "'\n", sep = "")
  
  sample_ids <- df[[1]]
  feature_matrix <- df[, -c(1, condition_col_index), drop = FALSE]
  
  # Transpose: rows = features/peaks, columns = samples
  transposed <- as.data.frame(t(feature_matrix), check.names = FALSE)
  colnames(transposed) <- sample_ids
  transposed <- cbind(ID = rownames(transposed), transposed)
  rownames(transposed) <- NULL
  
  cat("[PRAMIGO prep] Transposed shape: ", nrow(transposed), " rows x ", ncol(transposed),
      " columns (features x [ID + samples])\n", sep = "")
  cat("[PRAMIGO prep] First 3 row IDs (features): ",
      paste(head(transposed$ID, 3), collapse = ", "), "\n", sep = "")
  cat("[PRAMIGO prep] First 3 sample column names: ",
      paste(head(colnames(transposed)[-1], 3), collapse = ", "), "\n", sep = "")
  
  # New path: replaces the last part of the filename (e.g. "_MOFA") with "_PRAMIGO"
  # Metabolomics_TEST1_MOFA.tsv -> Metabolomics_TEST1_PRAMIGO.tsv
  new_path <- sub("_[^_]+\\.tsv$", "_PRAMIGO.tsv", original_path)
  
  write.table(transposed, new_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("[PRAMIGO prep] Written to: ", new_path, "\n\n", sep = "")
  
  new_path
}


# --- Resolve the correct Integration/INPUT file for PRAMIGO, mirroring the
# same path convention MOFA already uses, then apply the format fix
# (transpose + drop condition column) only where needed -------------------

get_pramigo_input_path <- function(i) {
  
  data_type <- COMMAND$DATA_TYPE[i]
  label     <- COMMAND$LABEL[i]
  
  if (data_type == "Transcriptomics") {
    
    integration_dir  <- paste(shared_dir, "/Integration/INPUT/", label, "_",
                              args[1], "_vs_", args[2], sep = "")
    integration_path <- paste(integration_dir, "/", label, "_",
                              args[1], "_vs_", args[2],
                              "_normalized_vst_variance.tsv", sep = "")
    
  } else if (data_type == "Methylomics") {
    
    integration_dir  <- paste(shared_dir, "/Integration/INPUT/Methylome_", label, "_",
                              args[1], "_vs_", args[2], sep = "")
    integration_path <- paste(integration_dir, "/", label, "_matrix_MOFA.tsv", sep = "")
    
  } else if (data_type == "Metabolomics") {
    
    integration_dir  <- paste(shared_dir, "/Integration/INPUT/Metabolomics_", label, "_",
                              args[1], "_vs_", args[2], sep = "")
    integration_path <- paste(integration_dir, "/Metabolomics_", label, "_MOFA.tsv", sep = "")
    
  } else if (data_type == "Undefined") {
    
    integration_dir  <- paste(shared_dir, "/Integration/INPUT/Undefined_", label, "_",
                              args[1], "_vs_", args[2], sep = "")
    integration_path <- paste(integration_dir, "/Undefined_", label, "_MOFA.tsv", sep = "")
    
  } else {
    stop("PRAMIGO: unknown DATA_TYPE '", data_type, "' - no Integration/INPUT convention defined.")
  }
  
  cat("[PRAMIGO prep] Resolved Integration/INPUT path for row ", i, " (", data_type, "): ",
      integration_path, "\n", sep = "")
  
  if (!file.exists(integration_path)) {
    stop("PRAMIGO: expected preprocessed file not found: ", integration_path)
  }
  
  # <-- QUESTO è l'UNICO pezzo cambiato rispetto a prima:
  if (data_type %in% c("Metabolomics", "Undefined")) {
    prepare_pramigo_matrix(integration_path)
  } else if (data_type %in% c("Transcriptomics", "Methylomics")) {
    prepare_pramigo_variance_filtered_matrix(integration_path)
  } else {
    integration_path
  }
}




prepare_pramigo_variance_filtered_matrix <- function(original_path, top_n = 5000) {
  
  cat("\n[PRAMIGO prep] Reading variance-annotated matrix: ", original_path, "\n", sep = "")
  
  df <- read.delim(original_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  cat("[PRAMIGO prep] Original shape: ", nrow(df), " rows x ", ncol(df),
      " columns ([ID + samples + variance])\n", sep = "")
  
  variance_col <- ncol(df)
  cat("[PRAMIGO prep] Last column ('", colnames(df)[variance_col],
      "') treated as the variance score\n", sep = "")
  
  df_sorted <- df[order(-df[[variance_col]]), , drop = FALSE]
  df_top    <- head(df_sorted, top_n)
  df_top    <- df_top[, -variance_col, drop = FALSE]
  
  cat("[PRAMIGO prep] Filtered to top ", top_n, " most variable rows -> shape: ",
      nrow(df_top), " rows x ", ncol(df_top), " columns ([ID + samples])\n", sep = "")
  cat("[PRAMIGO prep] First 3 IDs kept: ",
      paste(head(df_top[[1]], 3), collapse = ", "), "\n", sep = "")
  
  new_path <- sub("_[^_]+\\.tsv$", "_PRAMIGO.tsv", original_path)
  
  write.table(df_top, new_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("[PRAMIGO prep] Written to: ", new_path, "\n\n", sep = "")
  
  new_path
}




#if(args[1] == ""){
# taking input with showing the message
args = commandArgs(trailingOnly=TRUE)
#}

print(args)

#directory <- "/home/cristia/Scrivania/PhD/Bioinfo/Article_MULTI_BETA"
directory <- unlist(args[[3]])


# NEW: shared_dir is the mounted folder (e.g. "/shared" in Docker mode),
# passed as the 4th argument from app.R. Only the JSON lives there for now.
# Falls back to `directory` when not provided (e.g. local/non-Docker runs
# with only 3 args), so nothing breaks outside Docker mode.
shared_dir <- if (length(args) >= 4) unlist(args[[4]]) else directory


setwd(paste(directory,"_INSTALL",sep="/"))
renv::load(paste(directory,"_INSTALL",sep="/"))

setwd(directory)
print(directory)

library(vroom)


site_lib <- "/usr/local/lib/R/site-library"
if (dir.exists(site_lib) && !(site_lib %in% .libPaths())) {
  .libPaths(c(.libPaths(), site_lib))
}


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
    Cell_type <- as.character(COMMAND$LABEL[i])
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
    Cell_type <- as.character(COMMAND$LABEL[i])
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
    Cell_type <- as.character(COMMAND$LABEL[i])
    cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
    iterations = iterations + 1
    selection_samples = COMMAND$SELECTION[i]
    #purity_filter =COMMAND$PURITY[i]
    directory2 <- paste(directory,"/Transcriptomics/INPUT",sep="")
    #files <- grep("MATRIX*",list.files(directory2),value=TRUE)
    #files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
    
    
    if (!is_in_docker) {
      # --- Caso locale: comportamento originale, invariato ---
      source(paste(directory,"/Transcriptomics/Biomix_DGE_GENES_LIMMA.r",sep=""))
    } else {
      # --- Caso Docker: delega al container fratello ---
      runner_transcriptomics <- generateRunnerFunction(
        "ghcr.io/biomix-consortium/biomix-transcriptomics:latest",
        expose_qc_port = TRUE   # NEW — serve la porta 3839 per il QC preview Shiny
      )
      
      runner_transcriptomics(
        command = "Rscript",
        args    = c(
          paste0(directory, "/Transcriptomics/Biomix_DGE_GENES_LIMMA.r"),
          directory,
          Cell_type,
          args[[1]],
          args[[2]],
          shared_dir,
          iterations
        )
      )
    }
    
    cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
    gc()
    
  }else if (COMMAND$ANALYSIS[i] == "NO" & COMMAND$INTEGRATION[i] == "YES") {
    STATISTICS <- "NO"
    Cell_type <- as.character(COMMAND$LABEL[i])
    cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
    iterations = iterations + 1
    selection_samples = COMMAND$SELECTION[i]
    #purity_filter =COMMAND$PURITY[i]
    directory2 <- paste(directory,"/Transcriptomics/INPUT",sep="")
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
    Cell_type <- as.character(COMMAND$LABEL[i])
    cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
    iterations = iterations + 1
    selection_samples = COMMAND$SELECTION[i] # TO ADD
    #purity_filter =COMMAND$PURITY[i] #TO ADD
    directory2 <- paste(directory,"/Metabolomics/INPUT",sep="")
    files <- grep("MATRIX*",list.files(directory2),value=TRUE)
    files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
    if (!is_in_docker) {
      # --- Caso locale: comportamento originale, invariato ---
      source(paste(directory,"/Metabolomics/BiomiX_DMA.r",sep=""))
    } else {
      # --- Caso Docker: delega al container fratello ---
      #runner_transcriptomics <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-transcriptomics:latest")
      runner_metabolomics <- generateRunnerFunction(
        "ghcr.io/biomix-consortium/biomix-metabolomics:latest",
        expose_qc_port = TRUE   # NEW — serve la porta 3839 per il QC preview Shiny
      )
      
      
      runner_metabolomics(
        command = "Rscript",
        args    = c(
          paste0(directory, "/Metabolomics/BiomiX_DMA.r"),
          directory,
          Cell_type,
          args[[1]],
          args[[2]],
          shared_dir,        # NEW — così lo script metabolomics sa dove leggere il JSON
          iterations        # NEW — posizione 6, valore già corretto in questo punto del loop
          
          
        )
      )
    }
    #source(paste(directory,"/Metabolomics/BiomiX_DMA.r",sep = ""))
    cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
    gc()
    
  }else if (COMMAND$ANALYSIS[i] == "NO" & COMMAND$INTEGRATION[i] == "YES") {
    STATISTICS <- "NO"
    Cell_type <- as.character(COMMAND$LABEL[i])
    cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
    iterations = iterations + 1
    selection_samples = COMMAND$SELECTION[i] # TO ADD
    #purity_filter =COMMAND$PURITY[i] #TO ADD
    directory2 <- paste(directory,"/Metabolomics/INPUT",sep="")
    files <- grep("MATRIX*",list.files(directory2),value=TRUE)
    files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
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
    Cell_type <- as.character(COMMAND$LABEL[i])
    cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
    iterations = iterations + 1
    selection_samples = COMMAND$SELECTION[i] # TO ADD
    #purity_filter =COMMAND$PURITY[i] #TO ADD
    directory2 <- paste(directory,"/Methylomics/INPUT",sep="")
    files <- grep("MATRIX*",list.files(directory2),value=TRUE)
    files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
    if (!is_in_docker) {
      # --- Caso locale: comportamento originale, invariato ---
      source(paste(directory,"/Methylomics/BiomiX_DMA.r",sep=""))
    } else {
      # --- Caso Docker: delega al container fratello ---
      #runner_transcriptomics <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-transcriptomics:latest")
      runner_methylomics <- generateRunnerFunction(
        "ghcr.io/biomix-consortium/biomix-methylomics:latest",
        expose_qc_port = TRUE   # NEW — serve la porta 3839 per il QC preview Shiny
      )
      
      
      runner_methylomics(
        command = "Rscript",
        args    = c(
          paste0(directory, "/Methylomics/BiomiX_DMA.r"),
          directory,
          Cell_type,
          args[[1]],
          args[[2]],
          shared_dir,        # NEW — così lo script metabolomics sa dove leggere il JSON
          iterations        # NEW — posizione 6, valore già corretto in questo punto del loop
          
          
        )
      )
    }
    #source(paste(directory,"/Methylomics/BiomiX_DMA.r",sep=""))
    cat("\n\n\n\n\n  ", Cell_type, " analysis complete ^-^")
    gc()
    
  } else if (COMMAND$ANALYSIS[i] == "NO" & COMMAND$INTEGRATION[i] == "YES") {
    STATISTICS <- "NO"
    Cell_type <- as.character(COMMAND$LABEL[i])
    cat(paste( "\n\n\n\n\nStarting the ", Cell_type, " analysis \n\n\n\n\n", sep =""))
    iterations = iterations + 1
    selection_samples = COMMAND$SELECTION[i] # TO ADD
    #purity_filter =COMMAND$PURITY[i] #TO ADD
    directory2 <- paste(directory,"/Methylomics/INPUT",sep="")
    files <- grep("MATRIX*",list.files(directory2),value=TRUE)
    files_meta <- grep("METADATA*",list.files(directory2),value=TRUE)
    cat("\n\n\n\n\n  Reminder: \n The Single Omics Analysis was not selected for", Cell_type, "\n If you are performing multiomics analysis and you did not run it before please do it,\n without the multiomics pipeline will not find the input data")
    gc()
  }
}

#================================================================================

old_comp <-  Sys.time() # calculate time computation

if(COMMAND_MOFA[2,2] == "YES") {
  
  if(COMMAND_MOFA[1,2] == "MOFA_INTEGRATION") {        
    
    
      cat(paste( "\n\n\n\n\nStarting the MOFA analysis \n\n\n\n\n", sep =""))
      Sys.sleep(5)
      
      if (!is_in_docker) {
        # --- Caso locale: comportamento originale, invariato ---
        source(paste(directory,"/Integration/MOFA_MULTI2_AUTO.R",sep=""))
      } else {
        # --- Caso Docker: delega al container fratello ---
        #runner_transcriptomics <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-transcriptomics:latest")
        runner_MOFA <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-mofa-diablo:latest")
        
        runner_MOFA(
          command = "Rscript",
          args    = c(
            paste0(directory, "/Integration/MOFA_MULTI2_AUTO.R"),
            directory,
            args[[1]],
            args[[2]],
            shared_dir,        # NEW — così lo script metabolomics sa dove leggere il JSON
            iterations        # NEW — posizione 6, valore già corretto in questo punto del loop
            
            
          )
        )
        
      }

      cat("\n\n\n\n\n MOFA analysis complete ^-^")
    
  }
  
  if(COMMAND_MOFA[1,2] == "DIABLO_INTEGRATION") {
    
    cat(paste( "\n\n\n\n\nStarting the DIABLO analysis \n\n\n\n\n", sep =""))
    Sys.sleep(5)
    if (!is_in_docker) {
      # --- Caso locale: comportamento originale, invariato ---
      source(paste(directory,"/Integration/DIABLO_MULTI2.R",sep=""))
    } else {
      # --- Caso Docker: delega al container fratello ---
      #runner_transcriptomics <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-transcriptomics:latest")
      runner_DIABLO <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-mofa-diablo:latest")
      
      runner_DIABLO(
        command = "Rscript",
        args    = c(
          paste0(directory, "/Integration/DIABLO_MULTI2.R"),
          directory,
          Cell_type,
          args[[1]],
          args[[2]],
          shared_dir,        # NEW — così lo script metabolomics sa dove leggere il JSON
          iterations        # NEW — posizione 6, valore già corretto in questo punto del loop
          
          
        )
      )
    }
    #source(paste(directory,"/integration/DIABLO_MULTI2.R",sep=""))
    cat("\n\n\n\n\n DIABLO analysis complete ^-^")
    
  }
  
  if(COMMAND_MOFA[1,2] == "SNF_INTEGRATION") {
    
    cat(paste( "\n\n\n\n\nStarting the SNF analysis \n\n\n\n\n", sep =""))
    Sys.sleep(5)
    if (!is_in_docker) {
      # --- Caso locale: comportamento originale, invariato ---
      source(paste(directory,"/Integration/SNF_int.R",sep=""))
    } else {
      # --- Caso Docker: delega al container fratello ---
      #runner_transcriptomics <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-transcriptomics:latest")
      runner_SNF <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-snf-nemo:latest")
      
      runner_SNF(
        command = "Rscript",
        args    = c(
          paste0(directory, "/Integration/SNF_int.R"),
          directory,
          Cell_type,
          args[[1]],
          args[[2]],
          shared_dir,        # NEW — così lo script metabolomics sa dove leggere il JSON
          iterations        # NEW — posizione 6, valore già corretto in questo punto del loop
          
          
        )
      )
    }
    #source(paste(directory,"/integration/SNF_int.R",sep=""))
    cat("\n\n\n\n\n SNF analysis complete ^-^")
    
  }
  
  if(COMMAND_MOFA[1,2] == "NEMO_INTEGRATION") {
    
    cat(paste( "\n\n\n\n\nStarting the NEMO analysis \n\n\n\n\n", sep =""))
    Sys.sleep(5)
    
    if (!is_in_docker) {
      # --- Caso locale: comportamento originale, invariato ---
      source(paste(directory,"/Integration/NEMO_int.R",sep=""))
    } else {
      # --- Caso Docker: delega al container fratello ---
      #runner_transcriptomics <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-transcriptomics:latest")
      runner_NEMO <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-snf-nemo:latest")
      
      runner_NEMO(
        command = "Rscript",
        args    = c(
          paste0(directory, "/Integration/NEMO_int.R"),
          directory,
          Cell_type,
          args[[1]],
          args[[2]],
          shared_dir,        # NEW — così lo script metabolomics sa dove leggere il JSON
          iterations        # NEW — posizione 6, valore già corretto in questo punto del loop
          
          
        )
      )
    }
    
    #source(paste(directory,"/integration/NEMO_int.R",sep=""))
    cat("\n\n\n\n\n NEMO analysis complete ^-^")
    
  }
  
  
  if(COMMAND_MOFA[1,2] == "PRAMIGO_INTEGRATION") {
    
    
    cat("[PRAMIGO debug] is_in_docker =", is_in_docker, "\n")
    cat("[PRAMIGO debug] shared_dir   =", shared_dir, "\n")
    cat("[PRAMIGO debug] directory    =", directory, "\n")
    
    
    cat("\n\n\n\n\nStarting the PRAMIGO analysis \n\n\n\n\n")
    Sys.sleep(5)
    
    if (!is_in_docker) {
      # --- Caso locale: da definire quando servirà ---
    } else {
      
      idx_integration <- which(COMMAND$INTEGRATION == "YES")
      
      if (length(idx_integration) < 2) {
        cat("\n\n\n\n\n PRAMIGO requires at least two dataset selected for the integration — skipped.\n\n\n\n\n")
      } else {
        
        result_dir <- paste0(
          shared_dir, "/Integration/OUTPUT/PRAMIGO_",
          paste(COMMAND$DATA_TYPE[idx_integration], collapse = "_vs_")
        )
        dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
        
        # --- Estrazione parametri avanzati PRAMIGO dal JSON ----------------------
        # Mappatura confermata (3 righe x 6 colonne ADVANCED_OPTION_PRAMIGO_1.._6):
        #   riga1: epoch | n_batch | lr | resume | leiden_resolution | checkpoint_every
        #   riga2: n_hid | ae_n_hid | n_heads | n_layers | leiden_min_cluster_size | (eco resume)
        #   riga3: topology | omic1_corr_cutoff | omic2_corr_cutoff | sample_mode | cuda | (vuoto)
        adv <- COMMAND_ADVANCED
        
        PRAMIGO_EPOCH      <- as.character(adv$ADVANCED_OPTION_PRAMIGO_1[1])
        PRAMIGO_N_BATCH    <- as.character(adv$ADVANCED_OPTION_PRAMIGO_2[1])
        PRAMIGO_LR         <- as.character(adv$ADVANCED_OPTION_PRAMIGO_3[1])
        PRAMIGO_RESUME     <- as.character(adv$ADVANCED_OPTION_PRAMIGO_4[1])
        PRAMIGO_LEIDEN_RES <- as.character(adv$ADVANCED_OPTION_PRAMIGO_5[1])
        PRAMIGO_CHECKPOINT <- as.character(adv$ADVANCED_OPTION_PRAMIGO_6[1])
        
        PRAMIGO_N_HID      <- as.character(adv$ADVANCED_OPTION_PRAMIGO_1[2])
        PRAMIGO_AE_N_HID   <- as.character(adv$ADVANCED_OPTION_PRAMIGO_2[2])
        PRAMIGO_N_HEADS    <- as.character(adv$ADVANCED_OPTION_PRAMIGO_3[2])
        PRAMIGO_N_LAYERS   <- as.character(adv$ADVANCED_OPTION_PRAMIGO_4[2])
        PRAMIGO_LEIDEN_MIN <- as.character(adv$ADVANCED_OPTION_PRAMIGO_5[2])
        
        PRAMIGO_TOPOLOGY     <- as.character(adv$ADVANCED_OPTION_PRAMIGO_1[3])
        PRAMIGO_OMIC1_CUTOFF <- as.character(adv$ADVANCED_OPTION_PRAMIGO_2[3])
        PRAMIGO_OMIC2_CUTOFF <- as.character(adv$ADVANCED_OPTION_PRAMIGO_3[3])
        PRAMIGO_SAMPLE_MODE  <- as.character(adv$ADVANCED_OPTION_PRAMIGO_4[3])
        PRAMIGO_CUDA_PREF    <- as.character(adv$ADVANCED_OPTION_PRAMIGO_5[3])  # scelta esplicita utente
        
        # --cuda: 0 = usa GPU0, qualunque altro valore = CPU (confermato da
        # src/biomix_hgt.py: parser.add_argument('--cuda', type=int, default=1)).
        # Nessuna auto-detection lato R: si usa direttamente la scelta
        # dell'utente in UI. Il fallback a CPU se la GPU risulta assente resta
        # comunque gestito dentro biomix_hgt.py stesso (vedi il warning
        # "Falling back to CPU" osservato nei log).
        requested_cuda <- suppressWarnings(as.integer(PRAMIGO_CUDA_PREF))
        if (is.na(requested_cuda)) requested_cuda <- 1L  # default sicuro: CPU
        use_gpu <- (requested_cuda == 0L)
        
        cat(if (use_gpu) "\n\n GPU requested by the user - PRAMIGO will use CUDA.\n\n"
            else "\n\n CPU requested by the user - PRAMIGO will use the CPU.\n\n")
        
        runner_PRAMIGO <- generateRunnerFunction("biomix-pramigo:test", gpu = use_gpu) #ghcr.io/biomix-consortium/biomix-pramigo:latest
        
        cuda_arg <- c("--cuda", as.character(requested_cuda))
        
        # --- Parametri comuni a entrambi i casi (2 omiche / manifest) -----------
        common_args <- c(
          "--epoch", PRAMIGO_EPOCH,
          "--n_batch", PRAMIGO_N_BATCH,
          "--lr", PRAMIGO_LR,
          "--leiden_resolution", PRAMIGO_LEIDEN_RES,
          "--checkpoint_every", PRAMIGO_CHECKPOINT,
          "--n_hid", PRAMIGO_N_HID,
          "--ae_n_hid", PRAMIGO_AE_N_HID,
          "--n_heads", PRAMIGO_N_HEADS,
          "--n_layers", PRAMIGO_N_LAYERS,
          "--leiden_min_cluster_size", PRAMIGO_LEIDEN_MIN,
          "--topology", PRAMIGO_TOPOLOGY,
          "--sample_mode", PRAMIGO_SAMPLE_MODE
        )
        
        if (!is.na(PRAMIGO_RESUME) && !(PRAMIGO_RESUME %in% c("off", ""))) {
          common_args <- c(common_args, "--resume", PRAMIGO_RESUME)
        }
        
        if (length(idx_integration) == 2) {
          
          # --- Caso semplice: due omiche, CLI diretta ---
          pramigo_args <- c(
            "-u",
            "src/biomix_hgt.py",
            "--metadata_path", DIR_METADATA,
            "--omic1_path", get_pramigo_input_path(idx_integration[1]),
            "--omic2_path", get_pramigo_input_path(idx_integration[2]),
            "--result_dir", result_dir,
            "--omic1_corr_cutoff", PRAMIGO_OMIC1_CUTOFF,
            "--omic2_corr_cutoff", PRAMIGO_OMIC2_CUTOFF,
            common_args,
            cuda_arg
          )
          
        } else {
          
          # --- Caso generale: 3+ omiche, manifest.yaml scritto in shared_dir ---
          manifest_path <- paste0(result_dir, "/manifest.yaml")
          
          uniform_corr_cutoff <- as.numeric(PRAMIGO_OMIC1_CUTOFF)
          
          omics_list <- lapply(idx_integration, function(i) {
            list(
              name        = tolower(COMMAND$DATA_TYPE[i]),
              path        = get_pramigo_input_path(i),
              corr_cutoff = uniform_corr_cutoff
            )
          })
          
          manifest <- list(
            topology    = PRAMIGO_TOPOLOGY,
            sample_mode = PRAMIGO_SAMPLE_MODE,
            omics       = omics_list
          )
          
          yaml::write_yaml(manifest, manifest_path)
          
          pramigo_args <- c(
            "-u",
            "src/biomix_hgt.py",
            "--metadata_path", DIR_METADATA,
            "--omics_manifest", manifest_path,
            "--result_dir", result_dir,
            common_args,
            cuda_arg
          )
        }
        
        runner_PRAMIGO(
          command = "python",
          args    = pramigo_args
        )
      }
    }
    
    cat("\n\n\n\n\n PRAMIGO analysis complete ^-^")
  }
  
  
  
  
}


time_comput <-  Sys.time() - old_comp # calculate difference

#================================================================================

#Sys.sleep(5)

if(COMMAND_MOFA[1,2] == "DIABLO_INTEGRATION" | COMMAND_MOFA[1,2] == "MOFA_INTEGRATION" & COMMAND_MOFA[2,2] == "YES") {
  
  cat(paste( "\n\n\n\n\nStarting the integrated factors correlation vs clinical factors \n\n\n\n\n", sep =""))
  
  
  
  
  if (!is_in_docker) {
    # --- Caso locale: comportamento originale, invariato ---
    source(paste(directory,"/Clinical_data/BiomiX_Clinical.R",sep=""))
  } else {
    # --- Caso Docker: delega al container fratello ---
    #runner_transcriptomics <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-transcriptomics:latest")
    runner_Interpretation <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-interpretation:latest")
    
    runner_Interpretation(
      command = "Rscript",
      args    = c(
        paste0(directory, "/Clinical_data/BiomiX_Clinical.R"),
        directory,
        Cell_type,
        args[[1]],
        args[[2]],
        shared_dir,        # NEW — così lo script metabolomics sa dove leggere il JSON
        iterations        # NEW — posizione 6, valore già corretto in questo punto del loop
        
        
      )
    )
  }
  
  #source(paste(directory,"/Clinical_data/BiomiX_Clinical.R",sep=""))
  cat("\n\n\n\n\n integrated factors correlation vs clinical factors completed ^-^")
}



#================================================================================

#Sys.sleep(5)

if(COMMAND_MOFA[1,2] == "DIABLO_INTEGRATION" | COMMAND_MOFA[1,2] == "MOFA_INTEGRATION" & COMMAND_MOFA[2,2] == "YES") {
  
  cat(paste( "\n\n\n\n\nStarting the integrated factors pathway analysis \n\n\n\n\n", sep =""))
  
  if (!is_in_docker) {
    source(paste(directory,"/Integration/BiomiX_MOFA_MINER.R",sep=""))
  } else {
    runner_MOFA_MINER <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-interpretation:latest")
    
    runner_MOFA_MINER(
      command = "Rscript",
      args    = c(
        paste0(directory, "/Integration/BiomiX_MOFA_MINER.R"),
        directory,
        Cell_type,
        args[[1]],
        args[[2]],
        shared_dir,
        iterations
      )
    )
  }
  
  cat("\n\n\n\n\n integrated factors pathway analysis completed ^-^")
}


#================================================================================



#Sys.sleep(5)

if(COMMAND_MOFA[1,2] == "DIABLO_INTEGRATION" | COMMAND_MOFA[1,2] == "MOFA_INTEGRATION" & COMMAND_MOFA[2,2] == "YES") {
  
  
  cat(paste( "\n\n\n\n\nStarting the integrated factors articles matching \n\n\n\n\n", sep =""))
  
  if (!is_in_docker) {
    source(paste(directory,"/Integration/BiomiX_PUBMED.R",sep=""))
  } else {
    runner_PUBMED <- generateRunnerFunction("ghcr.io/biomix-consortium/biomix-interpretation:latest")
    
    runner_PUBMED(
      command = "Rscript",
      args    = c(
        paste0(directory, "/Integration/BiomiX_PUBMED.R"),
        directory,
        Cell_type,
        args[[1]],
        args[[2]],
        shared_dir,
        iterations
      )
    )
  }
  
  cat("\n\n\n\n\n integrated factors articles matching completed ^-^")
}


#================================================================================


#args[1] = "ICPLUS"
cat(paste( "\n\n\n\n\n Saving the files in the directory \n\n\n\n\n", sep =""))

Sys.sleep(5)

#directory2 <- paste(shared_dir, "/",COMMAND$DATA_TYPE[1], "/OUTPUT/",COMMAND$LABEL[1], "_",args[1], "_vs_", args[2], sep ="")
#print(directory2)
#setwd(directory2)

#COPY OMICS DIRECTORY

if (shared_dir != DIR_METADATA_output){
  
  for(output in 1:length(na.omit(COMMAND$LABEL))) {
    print(output)
    print(paste(shared_dir,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))
    if(file.exists(paste(shared_dir,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))){
      dir.create(path =  paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      file.copy(paste(shared_dir,"/",COMMAND$DATA_TYPE[output], "/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""), paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
      print(paste(DIR_METADATA_output,"/",COMMAND$DATA_TYPE[output],"/OUTPUT/",COMMAND$LABEL[output], "_",args[1], "_vs_", args[2], sep =""))
    }
  }
  
  if(COMMAND_MOFA[2,2] == "YES") {
    
    #COPY MOFA FILES
    output<-list.files(path = paste(shared_dir,"/","MOFA", "/OUTPUT/",sep =""), pattern = paste("^","AUTOMATIC_SEARCH_BEST_FACTORS_",args[1], "_vs_", args[2], sep=""))
    if(file.exists(paste(shared_dir,"/","MOFA", "/OUTPUT/",output, sep =""))){
      dir.create(path =  paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
      file.copy(paste(shared_dir,"/","MOFA", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep =""), overwrite = TRUE, recursive = FALSE, copy.mode = TRUE) 
      print(paste(output,"copied", sep = " "))
    }
    
    #COPY MOFA DIRECTORY
    file_to_transfer<-dir(path = paste(shared_dir,"/","MOFA", "/OUTPUT/",sep =""), pattern = paste("^","MOFA_",args[1], "_vs_", args[2], sep="") )
    for(output in file_to_transfer){
      if(file.exists(paste(shared_dir,"/","MOFA", "/OUTPUT/",output, sep =""))){
        dir.create(path =  paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        file.copy(paste(shared_dir,"/","MOFA", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","MOFA","/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
        #unlink(paste(shared_dir,"/","MOFA", "/OUTPUT/",output, sep =""), recursive=TRUE)
        print(paste(output,"copied", sep = " "))
      }}
    
    if(COMMAND_ADVANCED[1,12] == "YES" | COMMAND_ADVANCED[2,12]  == "YES") {
      
      output<-dir(path = paste(shared_dir,"/","Clinical_data", "/OUTPUT/",sep =""), pattern = paste("^","Clinical_",args[1], "_vs_", args[2], sep="") )
      if(file.exists(paste(shared_dir,"/","Clinical_data", "/OUTPUT/",output, sep =""))){
        dir.create(path =  paste(DIR_METADATA_output,"/","Clinical_data","/OUTPUT/", sep ="") ,  showWarnings = TRUE, recursive = TRUE, mode = "0777")
        file.copy(paste(shared_dir,"/","Clinical_data", "/OUTPUT/",output, sep =""), paste(DIR_METADATA_output,"/","Clinical_data","/OUTPUT/", sep =""), overwrite = TRUE, recursive = TRUE, copy.mode = TRUE) 
        #unlink(paste(shared_dir,"/","Clinical_data", "/OUTPUT/",output, sep =""), recursive=TRUE)
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



# print elapsed time
new <- Sys.time() - old # calculate difference

print(paste("Total time running:",new)) # print in nice format
print(paste("Total time integration:",time_comput)) # print in nice format

null_to_na <- function(x) {
  if (is.null(x) || length(x) == 0) NA_character_ else as.character(x)
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
    null_to_na(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2),
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2,
    COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3[1],
    null_to_na(COMMAND_ADVANCED$ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2),
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
    "NEMO_N°Neighbours", "NEMO_SNF_Scaling", "NEMO_only_common_samples",
    "SNF_NEMO_Enrichment_Variable", "SNF_NEMO_Survival_Variable", "SNF_NEMO_Ground_truth_variable",
    "SNF_NEMO_Max_cluster_expected", "SNF_NEMO_N°feature_heatmap", "SNF_NEMO_N°input_features",
    "DIABLO_Design_matrix_path", "Diablo_features_selection?", "Diablo_N°iterations"
    
  ),
  values = c(
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA,
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY_2[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_CLINIC_DATA_DIRECTORY[1:2],
    COMMAND_ADVANCED$ADVANCED_OPTION_SNF_OPTIONS,
    COMMAND_ADVANCED$ADVANCED_OPTION_NEMO_OPTIONS,
    COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES,
    COMMAND_ADVANCED$ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS,
    COMMAND_ADVANCED$ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN[1],
    COMMAND_ADVANCED$ADVANCED_OPTION_DIABLO_OPTIONS[1:2]
  )
)


matr_5 <- make_param_table(
  names = c(
    "Total_Time_running_seconds", "Time_integration_seconds"
    
  ),
  values = c(
    new,
    time_comput
  )
)



# ---------- OUTPUT DIRECTORY ----------

dir_out <- file.path(shared_dir, "Report_parameters")
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

addWorksheet(wb, "Running_time")
writeData(wb, "Running_time", matr_5)

saveWorkbook(wb, excel_file, overwrite = TRUE)


