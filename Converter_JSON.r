# Load required libraries
library(jsonlite)
library(readr)
library(dplyr)
library(tidyverse)

# Function to check mandatory fields
check_mandatory_fields <- function(data, fields, file_name) {
        missing_fields <- setdiff(fields, colnames(data))
        if (length(missing_fields) > 0) {
                stop(paste("Error: Missing mandatory fields in", file_name, "->", paste(missing_fields, collapse=", ")))
        }
}

### --- Process COMMANDS.tsv --- ###
commands <- read_tsv("COMMANDS.tsv", show_col_types = FALSE)
mandatory_fields_commands <- c("INDEX", "ANALYSIS", "DATA_TYPE", "INTEGRATION", 
                               "LABEL", "SELECTION", "DIRECTORIES", "PREVIEW")
check_mandatory_fields(commands, mandatory_fields_commands, "COMMANDS.tsv")
# commands <- rename(commands, DIRECTORIES = DIRECTORIES)
commands <- commands %>% select(2:ncol(commands))

### --- Process COMMANDS_MOFA.tsv --- ###
mofa <- read_tsv("COMMANDS_MOFA.tsv", col_names = FALSE, skip = 1, show_col_types = FALSE)
if (ncol(mofa) != 2) {
        stop("Error: COMMANDS_MOFA.tsv must have exactly two columns (key-value format).")
}
colnames(mofa) <- c("Parameter", "Value")
mofa_list <- setNames(as.list(mofa$Value), mofa$Parameter)

### --- Process COMMANDS_ADVANCED.tsv --- ###
advanced <- read_tsv("COMMANDS_ADVANCED.tsv", show_col_types = FALSE)

### --- Process command-line arguments --- ###
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Example assumption: args = c("SLE", "CTRL", "Directory_BiomiX")
if (length(args) < 3) {
        stop("Error: Expected at least 3 command-line arguments (e.g., GROUP_1, GROUP_2, BASE_DIR)")
}

# Assign argument names
arg_list <- list(
        GROUP_1 = args[1],
        GROUP_2 = args[2],
        BASE_DIR = args[3]
)

### --- Process directory.txt and directory_out.txt --- ###
metadata_dir <- tryCatch({
        readLines("directory.txt", warn = FALSE)[1]
}, error = function(e) {
        stop("Error reading directory.txt: ", e$message)
})

output_dir <- tryCatch({
        readLines("directory_out.txt", warn = FALSE)[1]
}, error = function(e) {
        stop("Error reading directory_out.txt: ", e$message)
})

# Create a list for the directory info
directory_info <- list(
        METADATA_DIR = metadata_dir,
        OUTPUT_DIR = output_dir
)


### --- Combine all into a single list --- ###
combined_json <- list(
        COMMANDS = commands,
        COMMANDS_MOFA = as_tibble(cbind(nms = names(as.tibble(mofa_list)), t(as.tibble(mofa_list)))),
        COMMANDS_ADVANCED = advanced,
        COMMAND_LINE_ARGS = arg_list,
        DIRECTORY_INFO = directory_info
)

# Write to a single JSON file
write( jsonlite::toJSON(combined_json, pretty = TRUE, auto_unbox = TRUE, na = "null"), "COMBINED_COMMANDS.json")

print("Combined JSON created successfully!")
