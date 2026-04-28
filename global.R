# =============================================================================
# global.R
# Loaded ONCE at app startup.
# Contains: libraries, global constants, and utility functions shared
# across all modules.
# =============================================================================

library(shiny)
library(jsonlite)
library(readxl)
library(readr)
library(dplyr)
library(base64enc)

# -----------------------------------------------------------------------------
# LOGO: load BiomiX_logo3.png and encode as base64.
# This approach works regardless of how Shiny serves static files,
# because the image is embedded directly in the HTML.
# The file is looked for in www/ first, then in the app root directory.
# -----------------------------------------------------------------------------

.logo_path <- local({
  candidates <- c(
    file.path("www", "BiomiX_logo3.png"),
    "BiomiX_logo3.png"
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) found[1] else NULL
})

BIOMIX_LOGO_TAG <- if (!is.null(.logo_path)) {
  # Encode as base64 and embed directly in the img src
  img_data <- base64enc::base64encode(.logo_path)
  tags$img(
    src    = paste0("data:image/png;base64,", img_data),
    height = "80px",
    style  = "object-fit: contain;"
  )
} else {
  # Fallback: try the www/ URL path (standard Shiny static serving)
  tags$img(
    src    = "BiomiX_logo3.png",
    height = "80px",
    style  = "object-fit: contain;"
  )
}

BIOMIX_VERSION <- "2.5"

# Available omics data types (used in the 6-input dropdowns)
OMICS_TYPES <- c(
  "Undefined",
  "Transcriptomics",
  "Proteomics",
  "Methylomics",
  "Metabolomics"
)

# Integration methods available
INTEGRATION_METHODS <- c("MOFA", "Diablo", "SNF", "NEMO", "MintTea")

# Index labels for the 6 inputs (used in the final JSON)
INPUT_INDICES <- paste0("input", 1:6)

# Default placeholder for empty/unused fields in the JSON
EMPTY_VALUE <- "X"


# -----------------------------------------------------------------------------
# COMMANDS_ADVANCED fixed rows 2 and 3
#
# The COMMANDS_ADVANCED array always contains exactly 3 rows:
#   Row 1 -> user-configured values (built dynamically in app.R)
#   Row 2 -> minimum/alternative thresholds (fixed, read by the pipeline)
#   Row 3 -> maximum values / extra options (fixed, read by the pipeline)
#
# Rows 2 and 3 are constants — they never change based on user input.
# They are defined here so the pipeline always finds the full 3-row structure.
# -----------------------------------------------------------------------------

COMMANDS_ADVANCED_ROW2 <- list(
  ADVANCED_OPTION_TRASCRIPTOMICS                      = "0.05",
  ADVANCED_OPTION_SUBGROUPING                         = "10",
  ADVANCED_OPTION_METABOLOMICS                        = 0.05,
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1         = "NO",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2       = NULL,
  ADVANCED_OPTION_METHYLOMICS                         = "0.05",
  ADVANCED_OPTION_MOFA                                = "fast",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2         = "10",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL     = "HMDB",
  ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY    = "300",
  ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY         = "10",
  ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL        = "NO",
  ADVANCED_OPTION_METADATA_FILTERING_1                = "Numerical",
  ADVANCED_OPTION_METADATA_FILTERING_2                = "Numerical",
  ADVANCED_OPTION_METADATA_FILTERING_3                = "Numerical",
  ADVANCED_OPTION_METADATA_FILTERING_4                = "Numerical",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX = "No",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES       = "X",
  ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY          = "X",
  ADVANCED_OPTION_CLINIC_DATA_DIRECTORY               = "X",
  ADVANCED_OPTION_CLUSTERING_OPTIONS                  = "ward.D2",
  ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY_2  = "10",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2       = NULL,
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3       = "X",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX = "No",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2   = "X",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4       = "15",
  ADVANCED_OPTION_SNF_OPTIONS                         = 0.5,
  ADVANCED_OPTION_NEMO_OPTIONS                        = "0.5",
  ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES          = "X",
  ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS            = 20,
  ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN             = "X",
  ADVANCED_OPTION_DIABLO_OPTIONS                      = "99"
)

COMMANDS_ADVANCED_ROW3 <- list(
  ADVANCED_OPTION_TRASCRIPTOMICS                      = "X",
  ADVANCED_OPTION_SUBGROUPING                         = "YES",
  ADVANCED_OPTION_METABOLOMICS                        = 3,
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1         = "15",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS1_2       = "hmdb/lipidmaps/metlin/kegg",
  ADVANCED_OPTION_METHYLOMICS                         = "450K",
  ADVANCED_OPTION_MOFA                                = "0.5",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2         = "1 (priority)/2 (priority)/3 (priority)",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_GENERAL     = "X",
  ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY    = "5000",
  ADVANCED_OPTION_MOFA_INTERPRETATION_PATHWAY         = "X",
  ADVANCED_OPTION_MOFA_INTERPRETATION_CLINICAL        = "X",
  ADVANCED_OPTION_METADATA_FILTERING_1                = NULL,
  ADVANCED_OPTION_METADATA_FILTERING_2                = NULL,
  ADVANCED_OPTION_METADATA_FILTERING_3                = NULL,
  ADVANCED_OPTION_METADATA_FILTERING_4                = NULL,
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_INDEX = "No",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES       = "X",
  ADVANCED_OPTION_METABOLOMICS_MS2_DIRECTORY          = "X",
  ADVANCED_OPTION_CLINIC_DATA_DIRECTORY               = "X",
  ADVANCED_OPTION_CLUSTERING_OPTIONS                  = "20",
  ADVANCED_OPTION_MOFA_INTERPRETATION_BIBLIOGRAPHY_2  = "X",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_2       = "rp",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3       = "X",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_3_INDEX = "No",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_FILES_MS2   = "X",
  ADVANCED_OPTION_METABOLOMICS_ANNOTATION_MS2_4       = "X",
  ADVANCED_OPTION_SNF_OPTIONS                         = 20,
  ADVANCED_OPTION_NEMO_OPTIONS                        = "X",
  ADVANCED_OPTION_SNF_NEMO_METADATA_FEATURES          = "CONDITION",
  ADVANCED_OPTION_SNF_NEMO_NUMERIC_OPTIONS            = 20,
  ADVANCED_OPTION_FILE_PATH_DIABLO_DESIGN             = "X",
  ADVANCED_OPTION_DIABLO_OPTIONS                      = "X"
)
# Supports .xlsx, .xls, .tsv, .csv, .txt
# Returns a dataframe or throws a user-readable error.
# -----------------------------------------------------------------------------

read_metadata <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  
  df <- tryCatch({
    if (ext %in% c("xlsx", "xls")) {
      read_excel(file_path)
    } else {
      read_delim(file_path, delim = "\t", show_col_types = FALSE)
    }
  }, error = function(e) {
    stop(paste("Error reading metadata file:", e$message))
  })
  
  # CONDITION column is mandatory
  if (!"CONDITION" %in% colnames(df)) {
    stop("The metadata file does not contain the required 'CONDITION' column.")
  }
  
  # Warn if CONDITION is numeric
  if (is.numeric(df$CONDITION)) {
    warning("The CONDITION column is numeric. Conditions should be character strings (e.g. 'Healthy', 'Disease').")
  }
  
  df
}


# -----------------------------------------------------------------------------
# UTILITY: extract unique conditions from metadata dataframe
# -----------------------------------------------------------------------------

get_conditions <- function(df) {
  as.character(unique(df$CONDITION))
}


# -----------------------------------------------------------------------------
# UTILITY: build the COMMANDS list (6 inputs) for the final JSON.
# Arguments:
#   inputs_data  = list collected from mod_inputs
# -----------------------------------------------------------------------------

build_commands <- function(inputs_data) {
  lapply(seq_len(6), function(i) {
    d <- inputs_data[[i]]
    list(
      INDEX       = INPUT_INDICES[i],
      ANALYSIS    = d$analysis,
      DATA_TYPE   = d$data_type,
      INTEGRATION = d$integration,
      LABEL       = if (nchar(trimws(d$label)) == 0) NULL else d$label,
      SELECTION   = d$selection,
      DIRECTORIES = if (nchar(trimws(d$path)) == 0) EMPTY_VALUE else d$path,
      PREVIEW     = d$preview
    )
  })
}


# -----------------------------------------------------------------------------
# UTILITY: build the COMMANDS_MOFA list for the final JSON.
# Arguments:
#   integration_data = list collected from mod_integration
# -----------------------------------------------------------------------------

build_commands_mofa <- function(integration_data) {
  list(
    list(nms = "Integration_type",              V2 = integration_data$integration_type),
    list(nms = "Data_Integration",              V2 = integration_data$data_integration),
    list(nms = "Number_of_factors_to_calculate",V2 = as.character(integration_data$n_factors)),
    list(nms = "Factor_to_explore",             V2 = as.character(integration_data$factor_explore)),
    list(nms = "Omics_overlay",                 V2 = as.character(integration_data$omics_overlay))
  )
}


# -----------------------------------------------------------------------------
# UTILITY: assemble and write the final combined JSON to disk.
# Arguments:
#   commands          = output of build_commands()
#   commands_mofa     = output of build_commands_mofa()
#   commands_advanced = list collected from all mod_adv_* modules
#   group1, group2    = condition / control strings
#   base_dir          = BiomiX base directory
#   metadata_dir      = path to the metadata file
#   output_dir        = results output directory
#   out_path          = destination path for the JSON file
# -----------------------------------------------------------------------------

write_combined_json <- function(commands,
                                commands_mofa,
                                commands_advanced,
                                group1,
                                group2,
                                base_dir,
                                metadata_dir,
                                output_dir,
                                out_path = "COMBINED_COMMANDS.json") {
  combined <- list(
    COMMANDS          = commands,
    COMMANDS_MOFA     = commands_mofa,
    COMMANDS_ADVANCED = commands_advanced,
    COMMAND_LINE_ARGS = list(
      GROUP_1  = group1,
      GROUP_2  = group2,
      BASE_DIR = base_dir
    ),
    DIRECTORY_INFO = list(
      METADATA_DIR = metadata_dir,
      OUTPUT_DIR   = output_dir
    )
  )
  
  write(
    jsonlite::toJSON(combined, pretty = TRUE, auto_unbox = TRUE, na = "null"),
    out_path
  )
  
  invisible(out_path)
}
