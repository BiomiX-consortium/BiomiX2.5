
print("looking for the BiomiX_folder")
print(getwd())

find_folder <- function(folder_name, depth = 5) {
        # Get the root directory
        root <- normalizePath("/")
        
        # Recursive function to search for the folder within the specified depth
        search_folder <- function(directory, current_depth) {
                if (current_depth <= depth) {
                        files <- list.files(directory, full.names = TRUE)
                        
                        matching_folders <- files[file.info(files)$isdir & basename(files) == folder_name]
                        
                        if (length(matching_folders) > 0) {
                                return(matching_folders[1])
                        } else {
                                subdirectories <- files[file.info(files)$isdir]
                                for (subdir in subdirectories) {
                                        result <- search_folder(subdir, current_depth + 1)
                                        if (!is.null(result)) {
                                                return(result)
                                        }
                                }
                        }
                }
                return(NULL)
        }
        
        # Call the recursive function starting from the root directory
        found_folder <- search_folder(root, 1)
        
        # Return the path to the matching folder (if found)
        if (!is.null(found_folder)) {
                return(found_folder)
        } else {
                return(NULL)
        }
}

setwd(find_folder("BiomiX2.5"))
setwd(paste(getwd(),"/", "_INSTALL", sep=""))

print("R Package checking..")
chooseCRANmirror(48, ind = TRUE)

print("R Package checking..")
print("getwd()")

if (file.exists("Package_linux.tar.gz") == TRUE){
        if (file.exists("Package_linux") == TRUE){
            print("R package already downloaded and decompressed")        
        }else{
            print("R package already downloaded, decompressing the .tar file ")   
            untar("Package_linux.tar.gz")
        }
}else{
        path<- paste(getwd(),"/", "Package_linux.tar.gz", sep="")
        options(timeout=6000)
        download.file("https://github.com/IxI-97/BiomiX/releases/download/v2.4/Package_linux.tar.gz", destfile = path, mode = "wb")
        untar("Package_linux.tar.gz")
}



setwd(paste(getwd(),"/", "Package_linux", sep=""))

chooseCRANmirror(48, ind = TRUE)

library('systemfonts')
options(Ncpus = 6)

install.packages("sys_3.4.2.tar.gz", repos=NULL, type="source")
install.packages("abind_1.4-5.tar.gz", repos=NULL, type="source")
install.packages("askpass_1.2.0.tar.gz", repos=NULL, type="source")
install.packages("backports_1.4.1.tar.gz", repos=NULL, type="source")
install.packages("base64enc_0.1-3.tar.gz", repos=NULL, type="source")
install.packages("BiocManager_1.30.24.tar.gz", repos=NULL, type="source")
install.packages("bit_4.0.5.tar.gz", repos=NULL, type="source")
install.packages("bit64_4.0.5.tar.gz", repos=NULL, type="source")
install.packages("glue_1.7.0.tar.gz", repos=NULL, type="source")
install.packages("cli_3.6.3.tar.gz", repos=NULL, type="source")
install.packages("rlang_1.1.4.tar.gz", repos=NULL, type="source")
install.packages("lifecycle_1.0.4.tar.gz", repos=NULL, type="source")
install.packages("vctrs_0.6.5.tar.gz", repos=NULL, type="source")
install.packages("blob_1.2.4.tar.gz", repos=NULL, type="source")
install.packages("boot_1.3-30.tar.gz", repos=NULL, type="source")
install.packages("brew_1.0-10.tar.gz", repos=NULL, type="source")
install.packages("fansi_1.0.6.tar.gz", repos=NULL, type="source")
install.packages("utf8_1.2.4.tar.gz", repos=NULL, type="source")
install.packages("pillar_1.9.0.tar.gz", repos=NULL, type="source")
install.packages("R6_2.5.1.tar.gz", repos=NULL, type="source")
install.packages("withr_3.0.1.tar.gz", repos=NULL, type="source")
install.packages("tidyselect_1.2.1.tar.gz", repos=NULL, type="source")

install.packages("magrittr_2.0.3.tar.gz", repos=NULL, type="source")
install.packages("pkgconfig_2.0.3.tar.gz", repos=NULL, type="source")
install.packages("tibble_3.2.1.tar.gz", repos=NULL, type="source")
install.packages("brio_1.1.5.tar.gz", repos=NULL, type="source")
install.packages("generics_0.1.3.tar.gz", repos=NULL, type="source")
install.packages("dplyr_1.1.4.tar.gz", repos=NULL, type="source")
install.packages("purrr_1.0.2.tar.gz", repos=NULL, type="source")
install.packages("stringi_1.8.4.tar.gz", repos=NULL, type="source")
install.packages("stringr_1.5.1.tar.gz", repos=NULL, type="source")
install.packages("cpp11_0.5.0.tar.gz", repos=NULL, type="source")
install.packages("tidyr_1.3.1.tar.gz", repos=NULL, type="source")

install.packages("hms_1.1.3.tar.gz", repos=NULL, type="source")
install.packages("clipr_0.8.0.tar.gz", repos=NULL, type="source")
install.packages("tzdb_0.4.0.tar.gz", repos=NULL, type="source")
install.packages("crayon_1.5.3.tar.gz", repos=NULL, type="source")
install.packages("prettyunits_1.2.0.tar.gz", repos=NULL, type="source")
install.packages("progress_1.2.3.tar.gz", repos=NULL, type="source")
install.packages("vroom_1.6.5.tar.gz", repos=NULL, type="source")
install.packages("readr_2.1.5.tar.gz", repos=NULL, type="source")
install.packages("forcats_1.0.0.tar.gz", repos=NULL, type="source")
install.packages("haven_2.5.4.tar.gz", repos=NULL, type="source")

install.packages("evaluate_0.24.0.tar.gz", repos=NULL, type="source")
install.packages("xfun_0.47.tar.gz", repos=NULL, type="source")
install.packages("highr_0.11.tar.gz", repos=NULL, type="source")
install.packages("yaml_2.3.10.tar.gz", repos=NULL, type="source")
install.packages("knitr_1.48.tar.gz", repos=NULL, type="source")

install.packages("fastmap_1.2.0.tar.gz", repos=NULL, type="source")
install.packages("cachem_1.1.0.tar.gz", repos=NULL, type="source")
install.packages("digest_0.6.37.tar.gz", repos=NULL, type="source")
install.packages("htmltools_0.5.8.1.tar.gz", repos=NULL, type="source")
install.packages("memoise_2.0.1.tar.gz", repos=NULL, type="source")
install.packages("mime_0.12.tar.gz", repos=NULL, type="source")
install.packages("fs_1.6.4.tar.gz", repos=NULL, type="source")
install.packages("rappdirs_0.3.3.tar.gz", repos=NULL, type="source")
install.packages("sass_0.4.9.tar.gz", repos=NULL, type="source")


install.packages("jquerylib_0.1.4.tar.gz", repos=NULL, type="source")
install.packages("jsonlite_1.8.8.tar.gz", repos=NULL, type="source")
install.packages("tinytex_0.52.tar.gz", repos=NULL, type="source")
install.packages("bslib_0.8.0.tar.gz", repos=NULL, type="source")
install.packages("fontawesome_0.5.2.tar.gz", repos=NULL, type="source")
install.packages("rmarkdown_2.28.tar.gz", repos=NULL, type="source")

install.packages("ps_1.7.7.tar.gz", repos=NULL, type="source")
install.packages("processx_3.8.4.tar.gz", repos=NULL, type="source")
install.packages("callr_3.7.6.tar.gz", repos=NULL, type="source")
install.packages("desc_1.4.3.tar.gz", repos=NULL, type="source")
install.packages("pkgbuild_1.4.4.tar.gz", repos=NULL, type="source")
install.packages("rprojroot_2.0.4.tar.gz", repos=NULL, type="source")
install.packages("pkgload_1.4.0.tar.gz", repos=NULL, type="source")
install.packages("praise_1.0.0.tar.gz", repos=NULL, type="source")
install.packages("diffobj_0.3.5.tar.gz", repos=NULL, type="source")
install.packages("rematch2_2.1.2.tar.gz", repos=NULL, type="source")
install.packages("waldo_0.5.2.tar.gz", repos=NULL, type="source")
install.packages("testthat_3.2.1.1.tar.gz", repos=NULL, type="source")


install.packages("labelled_2.13.0.tar.gz", repos=NULL, type="source")
install.packages("broom_1.0.6.tar.gz", repos=NULL, type="source")
install.packages("broom.helpers_1.16.0.tar.gz", repos=NULL, type="source")
install.packages("Deriv_4.1.3.tar.gz", repos=NULL, type="source")
install.packages("farver_2.1.2.tar.gz", repos=NULL, type="source")
install.packages("labeling_0.4.3.tar.gz", repos=NULL, type="source")
install.packages("colorspace_2.1-1.tar.gz", repos=NULL, type="source")
install.packages("munsell_0.5.1.tar.gz", repos=NULL, type="source")
install.packages("RColorBrewer_1.1-3.tar.gz", repos=NULL, type="source")
install.packages("viridisLite_0.4.2.tar.gz", repos=NULL, type="source")

install.packages("scales_1.3.0.tar.gz", repos=NULL, type="source")
install.packages("MASS_7.3-58.3.tar.gz", repos=NULL, type="source")
install.packages("isoband_0.2.7.tar.gz", repos=NULL, type="source")
install.packages("gtable_0.3.5.tar.gz", repos=NULL, type="source")
install.packages("lattice_0.22-6.tar.gz", repos=NULL, type="source")
install.packages("nlme_3.1-166.tar.gz", repos=NULL, type="source")
install.packages("Matrix_1.7-0.tar.gz", repos=NULL, type="source")
install.packages("mgcv_1.9-1.tar.gz", repos=NULL, type="source")
install.packages("ggplot2_3.5.1.tar.gz", repos=NULL, type="source")
install.packages("survival_3.7-0.tar.gz", repos=NULL, type="source")
install.packages("mvtnorm_1.2-6.tar.gz", repos=NULL, type="source")
install.packages("TH.data_1.1-2.tar.gz", repos=NULL, type="source")
install.packages("RUnit_0.4.33.tar.gz", repos=NULL, type="source")
install.packages("microbenchmark_1.4.10.tar.gz", repos=NULL, type="source")

install.packages("cowplot_1.1.3.tar.gz", repos=NULL, type="source")
install.packages("modelr_0.1.11.tar.gz", repos=NULL, type="source")
install.packages("doBy_4.6.22.tar.gz", repos=NULL, type="source")
install.packages("Rcpp_1.0.13.tar.gz", repos=NULL, type="source")
install.packages("minqa_1.2.8.tar.gz", repos=NULL, type="source")
install.packages("nloptr_2.1.1.tar.gz", repos=NULL, type="source")
install.packages("RcppEigen_0.3.4.0.1.tar.gz", repos=NULL, type="source")
install.packages("lme4_1.1-35.5.tar.gz", repos=NULL, type="source")
install.packages("numDeriv_2016.8-1.1.tar.gz", repos=NULL, type="source")
install.packages("pbkrtest_0.5.3.tar.gz", repos=NULL, type="source")

install.packages("carData_3.0-5.tar.gz", repos=NULL, type="source")
install.packages("nnet_7.3-19.tar.gz", repos=NULL, type="source")
install.packages("SparseM_1.84-2.tar.gz", repos=NULL, type="source")
install.packages("MatrixModels_0.5-3.tar.gz", repos=NULL, type="source")
install.packages("quantreg_5.98.tar.gz", repos=NULL, type="source")
install.packages("car_3.1-2.tar.gz", repos=NULL, type="source")
install.packages("plyr_1.8.9.tar.gz", repos=NULL, type="source")
install.packages("reshape2_1.4.4.tar.gz", repos=NULL, type="source")
install.packages("timeDate_4032.109.tar.gz", repos=NULL, type="source")
install.packages("codetools_0.2-20.tar.gz", repos=NULL, type="source")
install.packages("timechange_0.3.0.tar.gz", repos=NULL, type="source")
install.packages("lubridate_1.9.3.tar.gz", repos=NULL, type="source")
install.packages("clock_0.7.1.tar.gz", repos=NULL, type="source")
install.packages("gower_1.0.1.tar.gz", repos=NULL, type="source")

install.packages("rpart_4.1.23.tar.gz", repos=NULL, type="source")
install.packages("class_7.3-22.tar.gz", repos=NULL, type="source")
install.packages("data.table_1.15.4.tar.gz", repos=NULL, type="source")
install.packages("shape_1.4.6.1.tar.gz", repos=NULL, type="source")
install.packages("diagram_1.6.5.tar.gz", repos=NULL, type="source")
install.packages("KernSmooth_2.23-24.tar.gz", repos=NULL, type="source")
install.packages("globals_0.16.3.tar.gz", repos=NULL, type="source")
install.packages("listenv_0.9.1.tar.gz", repos=NULL, type="source")
install.packages("parallelly_1.38.0.tar.gz", repos=NULL, type="source")
install.packages("future_1.34.0.tar.gz", repos=NULL, type="source")
install.packages("future.apply_1.11.2.tar.gz", repos=NULL, type="source")
install.packages("progressr_0.14.0.tar.gz", repos=NULL, type="source")
install.packages("SQUAREM_2021.1.tar.gz", repos=NULL, type="source")
install.packages("lava_1.8.0.tar.gz", repos=NULL, type="source")


install.packages("prodlim_2024.06.25.tar.gz", repos=NULL, type="source")
install.packages("ipred_0.9-15.tar.gz", repos=NULL, type="source")
install.packages("hardhat_1.4.0.tar.gz", repos=NULL, type="source")
install.packages("recipes_1.1.0.tar.gz", repos=NULL, type="source")
install.packages("pROC_1.18.5.tar.gz", repos=NULL, type="source")
install.packages("ModelMetrics_1.2.2.2.tar.gz", repos=NULL, type="source")
install.packages("iterators_1.0.14.tar.gz", repos=NULL, type="source")
install.packages("foreach_1.5.2.tar.gz", repos=NULL, type="source")
install.packages("proxy_0.4-27.tar.gz", repos=NULL, type="source")
install.packages("e1071_1.7-14.tar.gz", repos=NULL, type="source")
install.packages("caret_6.0-94.tar.gz", repos=NULL, type="source")
install.packages("rematch_2.0.0.tar.gz", repos=NULL, type="source")
install.packages("cellranger_1.1.0.tar.gz", repos=NULL, type="source")
install.packages("GlobalOptions_0.1.2.tar.gz", repos=NULL, type="source")
install.packages("circlize_0.4.16.tar.gz", repos=NULL, type="source")
install.packages("commonmark_1.9.1.tar.gz", repos=NULL, type="source")
install.packages("conflicted_1.2.0.tar.gz", repos=NULL, type="source")
install.packages("corrplot_0.94.tar.gz", repos=NULL, type="source")
install.packages("openssl_2.2.1.tar.gz", repos=NULL, type="source")
install.packages("curl_5.2.1.tar.gz", repos=NULL, type="source")
install.packages("credentials_2.0.1.tar.gz", repos=NULL, type="source")
install.packages("DBI_1.2.3.tar.gz", repos=NULL, type="source")
install.packages("dbplyr_2.5.0.tar.gz", repos=NULL, type="source")


install.packages("rstudioapi_0.16.0.tar.gz", repos=NULL, type="source")
install.packages("zip_2.3.1.tar.gz", repos=NULL, type="source")
install.packages("gert_2.1.1.tar.gz", repos=NULL, type="source")
install.packages("gitcreds_0.1.2.tar.gz", repos=NULL, type="source")
install.packages("ini_0.3.1.tar.gz", repos=NULL, type="source")
install.packages("httr2_1.0.3.tar.gz", repos=NULL, type="source")
install.packages("gh_1.4.1.tar.gz", repos=NULL, type="source")
install.packages("whisker_0.4.1.tar.gz", repos=NULL, type="source")
install.packages("usethis_3.0.0.tar.gz", repos=NULL, type="source")
install.packages("ellipsis_0.3.2.tar.gz", repos=NULL, type="source")


install.packages("xtable_1.8-4.tar.gz", repos=NULL, type="source")
install.packages("sourcetools_0.1.7-1.tar.gz", repos=NULL, type="source")
install.packages("later_1.3.2.tar.gz", repos=NULL, type="source")
install.packages("promises_1.3.0.tar.gz", repos=NULL, type="source")
install.packages("httpuv_1.6.15.tar.gz", repos=NULL, type="source")
install.packages("shiny_1.9.1.tar.gz", repos=NULL, type="source")
install.packages("miniUI_0.1.1.1.tar.gz", repos=NULL, type="source")
install.packages("htmlwidgets_1.6.4.tar.gz", repos=NULL, type="source")
install.packages("profvis_0.3.8.tar.gz", repos=NULL, type="source")
install.packages("xopen_1.0.1.tar.gz", repos=NULL, type="source")
install.packages("sessioninfo_1.2.2.tar.gz", repos=NULL, type="source")
install.packages("rcmdcheck_1.4.0.tar.gz", repos=NULL, type="source")
install.packages("remotes_2.5.0.tar.gz", repos=NULL, type="source")
install.packages("xml2_1.3.6.tar.gz", repos=NULL, type="source")
install.packages("roxygen2_7.3.2.tar.gz", repos=NULL, type="source")
install.packages("rversions_2.1.2.tar.gz", repos=NULL, type="source")

install.packages("urlchecker_1.0.1.tar.gz", repos=NULL, type="source")
#install.packages("systemfonts_1.1.0.tar.gz", repos=NULL, type="source")
install.packages("textshaping_0.4.0.tar.gz", repos=NULL, type="source")


# IF A PACKAGE IS NOT WORKING ON DIFFERENT LINUX VERSION, YOU CAN ADD NEWER VERSION IN THE REPOSITORY TO INTALL THE IF THE STANDARD FAILS. 

pkg_files <- list.files(pattern = "ragg_", full.names = TRUE)

for (pkg in pkg_files) {
  message(paste("Trying:", pkg))
  install.packages(pkg, repos = NULL, type = "source")

  # Extract package name (assuming standard format like "mypackage_1.2.3.tar.gz")
  pkg_name <- sub("(_[0-9.]+)?\\.tar\\.gz$", "", basename(pkg))

  # Check if the package loads
  load_attempt <- try(suppressMessages(library(pkg_name, character.only = TRUE)), silent = TRUE)
  
  if (!inherits(load_attempt, "try-error")) {
    message(paste("Installed successfully:", pkg))
    break
  } else {
    message(paste("Installation failed:", pkg))
  }
}

#install.packages("ragg_1.3.2.tar.gz", repos=NULL, type="source")
#install.packages("ragg_1.3.3.tar.gz", repos=NULL, type="source")
install.packages("uuid_1.2-1.tar.gz", repos=NULL, type="source")
install.packages("downlit_0.4.4.tar.gz", repos=NULL, type="source")
install.packages("pkgdown_2.1.0.tar.gz", repos=NULL, type="source")
install.packages("devtools_2.4.5.tar.gz", repos=NULL, type="source")




install.packages("dtplyr_1.3.1.tar.gz", repos=NULL, type="source")
install.packages("httr_1.4.7.tar.gz", repos=NULL, type="source")
install.packages("WriteXLS_6.7.0.tar.gz", repos=NULL, type="source")
install.packages("rjson_0.2.21.tar.gz", repos=NULL, type="source")
install.packages("enrichR_3.2.tar.gz", repos=NULL, type="source")
install.packages("gargle_1.5.2.tar.gz", repos=NULL, type="source")
install.packages("broom.helpers_1.16.0.tar.gz", repos=NULL, type="source")
install.packages("patchwork_1.2.0.tar.gz", repos=NULL, type="source")
install.packages("ggstats_0.6.0.tar.gz", repos=NULL, type="source")
install.packages("GGally_2.2.1.tar.gz", repos=NULL, type="source")
install.packages("rstatix_0.7.2.tar.gz", repos=NULL, type="source")

install.packages("ggrepel_0.9.5.tar.gz", repos=NULL, type="source")
install.packages("ggsci_3.2.0.tar.gz", repos=NULL, type="source")
install.packages("ggsignif_0.6.4.tar.gz", repos=NULL, type="source")
install.packages("gridExtra_2.3.tar.gz", repos=NULL, type="source")
install.packages("polynom_1.4-1.tar.gz", repos=NULL, type="source")
install.packages("ggpubr_0.6.0.tar.gz", repos=NULL, type="source")


install.packages("googledrive_2.1.1.tar.gz", repos=NULL, type="source")
install.packages("ids_1.0.1.tar.gz", repos=NULL, type="source")
install.packages("googlesheets4_1.1.1.tar.gz", repos=NULL, type="source")
install.packages("here_1.0.1.tar.gz", repos=NULL, type="source")
#install.packages("igraph_1.4.2.zip", repos=NULL, type="source")
#install.packages("lattice_0.21-8.zip", repos=NULL, type="source")
#install.packages("Matrix_1.5-4.zip", repos=NULL, type="source")
install.packages("ncdf4_1.23.tar.gz", repos=NULL, type="source")

install.packages("png_0.1-8.tar.gz", repos=NULL, type="source")
install.packages("RcppTOML_0.2.2.tar.gz", repos=NULL, type="source")
install.packages("readxl_1.4.3.tar.gz", repos=NULL, type="source")
install.packages("XML_3.99-0.17.tar.gz", repos=NULL, type="source")
install.packages("rentrez_1.2.3.tar.gz", repos=NULL, type="source")
install.packages("reprex_2.1.1.tar.gz", repos=NULL, type="source")
install.packages("reshape_0.8.9.tar.gz", repos=NULL, type="source")
install.packages("reticulate_1.38.0.tar.gz", repos=NULL, type="source")
install.packages("rlist_0.4.6.2.tar.gz", repos=NULL, type="source")
install.packages("selectr_0.4-2.tar.gz", repos=NULL, type="source")
install.packages("rvest_1.0.4.tar.gz", repos=NULL, type="source")
install.packages("tidyverse_2.0.0.tar.gz", repos=NULL, type="source")


library("httr2")
library("devtools")
library("vroom")
library("lava")
library("recipes")
library("dplyr")
library("future.apply")
library("stringr")
library("circlize")
library("ggplot2")
library("ggrepel")
library("enrichR")
library("rlist")
library("tidyverse")
library("data.table")
library("caret")
library("reticulate")
library("XML")
library("xml2")
library("rentrez")
library("remotes")
library("GGally")
library("ggpubr")
library("htmltools")
library("ggpubr")
library("ncdf4")



install.packages("bitops_1.0-8.tar.gz", repos=NULL, type="source")
install.packages("formatR_1.14.tar.gz", repos=NULL, type="source")
install.packages("plogr_0.2.0.tar.gz", repos=NULL, type="source")
install.packages("RCurl_1.98-1.16.tar.gz", repos=NULL, type="source")
install.packages("zlibbioc_1.50.0.tar.gz", repos=NULL, type="source")
install.packages("matrixStats_1.3.0.tar.gz", repos=NULL, type="source")
install.packages("lambda.r_1.2.4.tar.gz", repos=NULL, type="source")
install.packages("futile.options_1.0.1.tar.gz", repos=NULL, type="source")
install.packages("RSQLite_2.3.7.tar.gz", repos=NULL, type="source")
install.packages("BiocGenerics_0.50.0.tar.gz", repos=NULL, type="source")
install.packages("S4Vectors_0.42.1.tar.gz", repos=NULL, type="source")
install.packages("IRanges_2.38.1.tar.gz", repos=NULL, type="source")

install.packages("GenomeInfoDbData_1.2.2.tar.gz", repos=NULL, type="source")
install.packages("UCSC.utils_1.0.0.tar.gz", repos=NULL, type="source")
install.packages("GenomeInfoDb_1.40.1.tar.gz", repos=NULL, type="source")
install.packages("XVector_0.44.0.tar.gz", repos=NULL, type="source")
install.packages("Biostrings_2.72.1.tar.gz", repos=NULL, type="source")
install.packages("KEGGREST_1.44.1.tar.gz", repos=NULL, type="source")
install.packages("MatrixGenerics_1.16.0.tar.gz", repos=NULL, type="source")
install.packages("futile.logger_1.4.3.tar.gz", repos=NULL, type="source")
install.packages("snow_0.4-4.tar.gz", repos=NULL, type="source")
install.packages("S4Arrays_1.4.1.tar.gz", repos=NULL, type="source")
install.packages("SparseArray_1.4.8.tar.gz", repos=NULL, type="source")
install.packages("DelayedArray_0.30.1.tar.gz", repos=NULL, type="source")
install.packages("Biobase_2.64.0.tar.gz", repos=NULL, type="source")
install.packages("AnnotationDbi_1.66.0.tar.gz", repos=NULL, type="source")
install.packages("annotate_1.82.0.tar.gz", repos=NULL, type="source")



install.packages("GenomicRanges_1.56.1.tar.gz", repos=NULL, type="source")
install.packages("BH_1.84.0-0.tar.gz", repos=NULL, type="source")
install.packages("BiocParallel_1.38.0.tar.gz", repos=NULL, type="source")
install.packages("genefilter_1.86.0.tar.gz", repos=NULL, type="source")
install.packages("locfit_1.5-9.10.tar.gz", repos=NULL, type="source")
install.packages("geneplotter_1.82.0.tar.gz", repos=NULL, type="source")
install.packages("RcppArmadillo_14.0.0-1.tar.gz", repos=NULL, type="source")
install.packages("DBI_1.2.3.tar.gz", repos=NULL, type="source")
install.packages("GenomeInfoDbData_1.2.12.tar.gz", repos=NULL, type="source")
install.packages("SummarizedExperiment_1.34.0.tar.gz", repos=NULL, type="source")
install.packages("BiocVersion_3.19.1.tar.gz", repos=NULL, type="source")
install.packages("DESeq2_1.44.0.tar.gz", repos=NULL, type="source")

#BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)
#BiocManager::install("DESeq2")
library(DESeq2)

install.packages("doParallel_1.0.17.tar.gz", repos=NULL, type="source")
install.packages("cluster_2.1.6.tar.gz", repos=NULL, type="source")
install.packages("clue_0.3-65.tar.gz", repos=NULL, type="source")
install.packages("GetoptLong_1.0.5.tar.gz", repos=NULL, type="source")
install.packages("ComplexHeatmap_2.20.0.tar.gz", repos=NULL, type="source")

library(ComplexHeatmap)
install.packages("Rhtslib_3.0.0.tar.gz", repos=NULL, type="source")
install.packages("Rsamtools_2.20.0.tar.gz", repos=NULL, type="source")
install.packages("GenomicAlignments_1.40.0.tar.gz", repos=NULL, type="source")
install.packages("BiocIO_1.14.0.tar.gz", repos=NULL, type="source")
install.packages("restfulr_0.0.15.tar.gz", repos=NULL, type="source")
install.packages("rtracklayer_1.64.0.tar.gz", repos=NULL, type="source")
install.packages("GenomicFeatures_1.56.0.tar.gz", repos=NULL, type="source")
install.packages("filelock_1.0.3.tar.gz", repos=NULL, type="source")
install.packages("BiocFileCache_2.12.0.tar.gz", repos=NULL, type="source")
install.packages("biomaRt_2.60.1.tar.gz", repos=NULL, type="source")
install.packages("TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2.tar.gz", repos=NULL, type="source")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

install.packages("rngtools_1.5.2.tar.gz", repos=NULL, type="source")
install.packages("R.methodsS3_1.8.2.tar.gz", repos=NULL, type="source")
install.packages("doRNG_1.8.6.tar.gz", repos=NULL, type="source")
install.packages("multtest_2.60.0.tar.gz", repos=NULL, type="source")
install.packages("scrime_1.3.5.tar.gz", repos=NULL, type="source")
install.packages("base64_2.0.1.tar.gz", repos=NULL, type="source")
install.packages("sparseMatrixStats_1.16.0.tar.gz", repos=NULL, type="source")
install.packages("Rhdf5lib_1.26.0.tar.gz", repos=NULL, type="source")
install.packages("rhdf5filters_1.16.0.tar.gz", repos=NULL, type="source")
install.packages("rhdf5_2.48.0.tar.gz", repos=NULL, type="source")
install.packages("beanplot_1.3.1.tar.gz", repos=NULL, type="source")
install.packages("statmod_1.5.0.tar.gz", repos=NULL, type="source")
install.packages("limma_3.60.4.tar.gz", repos=NULL, type="source")
install.packages("bumphunter_1.46.0.tar.gz", repos=NULL, type="source")
install.packages("siggenes_1.78.0.tar.gz", repos=NULL, type="source")
install.packages("preprocessCore_1.66.0.tar.gz", repos=NULL, type="source")
install.packages("illuminaio_0.46.0.tar.gz", repos=NULL, type="source")
install.packages("DelayedMatrixStats_1.26.0.tar.gz", repos=NULL, type="source")
install.packages("mclust_6.1.1.tar.gz", repos=NULL, type="source")
install.packages("quadprog_1.5-8.tar.gz", repos=NULL, type="source")
install.packages("R.oo_1.26.0.tar.gz", repos=NULL, type="source")
install.packages("R.utils_2.12.3.tar.gz", repos=NULL, type="source")
install.packages("GEOquery_2.72.0.tar.gz", repos=NULL, type="source")
install.packages("HDF5Array_1.32.1.tar.gz", repos=NULL, type="source")
install.packages("nor1mix_1.3-3.tar.gz", repos=NULL, type="source")
install.packages("minfi_1.50.0.tar.gz", repos=NULL, type="source")
install.packages("IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0.tar.gz", repos=NULL, type="source")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

install.packages("deldir_2.0-4.tar.gz", repos=NULL, type="source")
install.packages("lazyeval_0.2.2.tar.gz", repos=NULL, type="source")
install.packages("crosstalk_1.2.1.tar.gz", repos=NULL, type="source")
install.packages("DT_0.33.tar.gz", repos=NULL, type="source")
install.packages("interactiveDisplayBase_1.42.0.tar.gz", repos=NULL, type="source")
install.packages("AnnotationFilter_1.28.0.tar.gz", repos=NULL, type="source")
install.packages("ProtGenerics_1.36.0.tar.gz", repos=NULL, type="source")
install.packages("dichromat_2.0-0.1.tar.gz", repos=NULL, type="source")
install.packages("jpeg_0.1-10.tar.gz", repos=NULL, type="source")
install.packages("interp_1.1-6.tar.gz", repos=NULL, type="source")
install.packages("affyio_1.74.0.tar.gz", repos=NULL, type="source")
install.packages("BSgenome_1.72.0.tar.gz", repos=NULL, type="source")
install.packages("VariantAnnotation_1.50.0.tar.gz", repos=NULL, type="source")
install.packages("AnnotationHub_3.12.0.tar.gz", repos=NULL, type="source")
install.packages("gtools_3.9.5.tar.gz", repos=NULL, type="source")
install.packages("permute_0.9-7.tar.gz", repos=NULL, type="source")
install.packages("beachmat_2.20.0.tar.gz", repos=NULL, type="source")
install.packages("latticeExtra_0.6-30.tar.gz", repos=NULL, type="source")

install.packages("foreign_0.8-87.tar.gz", repos=NULL, type="source")
install.packages("Formula_1.2-5.tar.gz", repos=NULL, type="source")
install.packages("checkmate_2.3.2.tar.gz", repos=NULL, type="source")
install.packages("htmlTable_2.4.3.tar.gz", repos=NULL, type="source")
install.packages("viridis_0.6.5.tar.gz", repos=NULL, type="source")
install.packages("Hmisc_5.1-3.tar.gz", repos=NULL, type="source")
install.packages("ensembldb_2.28.0.tar.gz", repos=NULL, type="source")
install.packages("biovizBase_1.52.0.tar.gz", repos=NULL, type="source")

install.packages("nleqslv_3.3.5.tar.gz", repos=NULL, type="source")
install.packages("affy_1.82.0.tar.gz", repos=NULL, type="source")
install.packages("ExperimentHub_2.12.0.tar.gz", repos=NULL, type="source")
install.packages("gtools_3.9.5.tar.gz", repos=NULL, type="source")
install.packages("bsseq_1.40.0.tar.gz", repos=NULL, type="source")
install.packages("edgeR_4.2.1.tar.gz", repos=NULL, type="source")
install.packages("DSS_2.52.0.tar.gz", repos=NULL, type="source")
install.packages("Gviz_1.48.0.tar.gz", repos=NULL, type="source")
install.packages("ROC_1.80.0.tar.gz", repos=NULL, type="source")####
install.packages("BiasedUrn_2.0.12.tar.gz", repos=NULL, type="source")
install.packages("ruv_0.9.7.1.tar.gz", repos=NULL, type="source")

install.packages("org.Hs.eg.db_3.19.1.tar.gz", repos=NULL, type="source")
install.packages("FDb.InfiniumMethylation.hg19_2.2.0.tar.gz", repos=NULL, type="source")
install.packages("methylumi_2.50.0.tar.gz", repos=NULL, type="source")
install.packages("lumi_2.56.0.tar.gz", repos=NULL, type="source")
install.packages("fastICA_1.2-5.1.tar.gz", repos=NULL, type="source")
install.packages("JADE_2.0-4.tar.gz", repos=NULL, type="source")
install.packages("RPMM_1.25.tar.gz", repos=NULL, type="source")
install.packages("prettydoc_0.4.1.tar.gz", repos=NULL, type="source")

install.packages("globaltest_5.58.0.tar.gz", repos=NULL, type="source")
install.packages("sva_3.52.0.tar.gz", repos=NULL, type="source")
install.packages("DNAcopy_1.78.0.tar.gz", repos=NULL, type="source")
install.packages("impute_1.78.0.tar.gz", repos=NULL, type="source")
install.packages("marray_1.82.0.tar.gz", repos=NULL, type="source")
install.packages("txdbmaker_1.0.1.tar.gz", repos=NULL, type="source")
install.packages("geneLenDataBase_1.40.1.tar.gz", repos=NULL, type="source")
install.packages("GO.db_3.19.1.tar.gz", repos=NULL, type="source")
install.packages("goseq_1.56.0.tar.gz", repos=NULL, type="source")
install.packages("IlluminaHumanMethylation450kmanifest_0.4.0.tar.gz", repos=NULL, type="source")#
install.packages("IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz", repos=NULL, type="source")
install.packages("wateRmelon_2.10.0.tar.gz", repos=NULL, type="source")
install.packages("IlluminaHumanMethylationEPICmanifest_0.3.0.tar.gz", repos=NULL, type="source")
install.packages("missMethyl_1.38.0.tar.gz", repos=NULL, type="source")
install.packages("kpmt_0.1.0.tar.gz", repos=NULL, type="source")
install.packages("qvalue_2.36.0.tar.gz", repos=NULL, type="source")
install.packages("isva_1.9.tar.gz", repos=NULL, type="source")
install.packages("shinythemes_1.2.0.tar.gz", repos=NULL, type="source")
install.packages("dendextend_1.17.1.tar.gz", repos=NULL, type="source")
install.packages("combinat_0.0-8.tar.gz", repos=NULL, type="source")


install.packages("Illumina450ProbeVariants.db_1.40.0.tar.gz", repos=NULL, type="source")
install.packages("plotly_4.10.4.tar.gz", repos=NULL, type="source")
install.packages("ChAMPdata_2.36.0.tar.gz", repos=NULL, type="source")
install.packages("DMRcate_3.0.5.tar.gz", repos=NULL, type="source")
install.packages("ChAMP_2.34.0.tar.gz", repos=NULL, type="source")

library(ChAMP)

install.packages("sitmo_2.0.2.tar.gz", repos=NULL, type="source")
install.packages("FNN_1.1.4.tar.gz", repos=NULL, type="source")
install.packages("irlba_2.3.5.1.tar.gz", repos=NULL, type="source")
install.packages("RcppAnnoy_0.0.22.tar.gz", repos=NULL, type="source")
install.packages("RSpectra_0.16-2.tar.gz", repos=NULL, type="source")
install.packages("dqrng_0.4.1.tar.gz", repos=NULL, type="source")
install.packages("RcppProgress_0.4.2.tar.gz", repos=NULL, type="source")
install.packages("dir.expiry_1.12.0.tar.gz", repos=NULL, type="source")
install.packages("basilisk.utils_1.16.0.tar.gz", repos=NULL, type="source")
install.packages("pheatmap_1.0.12.tar.gz", repos=NULL, type="source")
install.packages("Rtsne_0.17.tar.gz", repos=NULL, type="source")
install.packages("uwot_0.2.2.tar.gz", repos=NULL, type="source")
install.packages("basilisk_1.16.0.tar.gz", repos=NULL, type="source")
install.packages("MOFA2_1.14.0.tar.gz", repos=NULL, type="source")

library(MOFA2)


install.packages("mzR_2.38.0.tar.gz", repos=NULL, type="source")
install.packages("MsCoreUtils_1.16.1.tar.gz", repos=NULL, type="source")
install.packages("vsn_3.72.0.tar.gz", repos=NULL, type="source")
install.packages("pcaMethods_1.96.0.tar.gz", repos=NULL, type="source")
install.packages("MALDIquant_1.22.3.tar.gz", repos=NULL, type="source")
install.packages("mzID_1.42.0.tar.gz", repos=NULL, type="source")
install.packages("BiocBaseUtils_1.6.0.tar.gz", repos=NULL, type="source")
install.packages("MultiAssayExperiment_1.30.3.tar.gz", repos=NULL, type="source")
install.packages("igraph_2.0.3.tar.gz ", repos=NULL, type="source")
install.packages("QFeatures_1.14.2.tar.gz", repos=NULL, type="source")

install.packages("PSMatch_1.8.0.tar.gz", repos=NULL, type="source")
install.packages("MSnbase_2.30.1.tar.gz", repos=NULL, type="source")
library(MSnbase)

install.packages("Rdisop_1.64.0.tar.gz", repos=NULL, type="source")

library(Rdisop)

library(devtools)
library(remotes)
library(igraph)



install.packages("RJSONIO_1.3-1.9.tar.gz", repos=NULL, type="source")
#devtools::install_github("lzyacht/cmmr", dependencies = FALSE)
install.packages("cmmr_1.0.5.tar.gz", repos=NULL, type="source")


library(cmmr)

install.packages("pbapply_1.7-2.tar.gz", repos=NULL, type="source")
#remotes::install_gitlab("jaspershen/masstools", dependencies = FALSE)
remotes::install_local("masstools_1.0.13.tar.gz", upgrade= "never")

library(masstools)

install.packages("openxlsx_4.2.6.1.tar.gz", repos=NULL, type="source")
#remotes::install_github("tidymass/massdataset", dependencies = FALSE)
remotes::install_local("massdataset_1.0.34.tar.gz", upgrade= "never")

library(massdataset)

install.packages("furrr_0.3.1.tar.gz", repos=NULL, type="source")
#remotes::install_github("tidymass/metid", dependencies = FALSE)
remotes::install_local("metid_1.2.35.tar.gz", upgrade= "never")

library(metid)

install.packages("polyclip_1.10-7.tar.gz", repos=NULL, type="source")
install.packages("tweenr_2.0.3.tar.gz", repos=NULL, type="source")
install.packages("graphlayouts_1.1.1.tar.gz", repos=NULL, type="source")
install.packages("tidygraph_1.3.1.tar.gz", repos=NULL, type="source")
install.packages("ggforce_0.4.2.tar.gz", repos=NULL, type="source")
install.packages("bookdown_0.40.tar.gz ", repos=NULL, type="source")
install.packages("BiocStyle_2.32.1.tar.gz", repos=NULL, type="source")
install.packages("ggraph_2.2.1.tar.gz", repos=NULL, type="source")
#remotes::install_gitlab("tidymass/metpath", dependencies = FALSE)
remotes::install_local("metpath_1.0.8.tar.gz", upgrade= "never")


library(metpath)



install.packages("slam_0.1-52.tar.gz", repos=NULL, type="source")
install.packages("NLP_0.3-0.tar.gz", repos=NULL, type="source")
install.packages("stringdist_0.9.12.tar.gz", repos=NULL, type="source")
install.packages("zoo_1.8-12.tar.gz", repos=NULL, type="source")
install.packages("tm_0.7-14.tar.gz", repos=NULL, type="source")
install.packages("ISOcodes_2024.02.12.tar.gz ", repos=NULL, type="source")
install.packages("stopwords_2.3.tar.gz", repos=NULL, type="source")
install.packages("synthesisr_0.3.0.tar.gz", repos=NULL, type="source")
install.packages("ngram_3.2.3.tar.gz", repos=NULL, type="source")
install.packages("SnowballC_0.7.1.tar.gz", repos=NULL, type="source")
install.packages("changepoint_2.2.4.tar.gz", repos=NULL, type="source")
#remotes::install_github("elizagrames/litsearchr", dependencies = FALSE)
remotes::install_local("litsearchr_1.0.0.tar.gz", upgrade= "never")

library(litsearchr)

install.packages("insight_0.20.3.tar.gz", repos=NULL, type="source")
install.packages("datawizard_0.12.2.tar.gz", repos=NULL, type="source")
install.packages("sjlabelled_1.2.0.tar.gz", repos=NULL, type="source")

#devtools::install_github("strengejacke/sjmisc", dependencies = FALSE)
remotes::install_local("sjmisc_2.8.10.tar.gz", upgrade= "never")

library(sjmisc)


install.packages("umap_0.2.10.0.tar.gz", repos=NULL, type="source")
install.packages("randomForest_4.7-1.2.tar.gz", repos=NULL, type="source")
install.packages("itertools_0.1-3.tar.gz", repos=NULL, type="source")
install.packages("missForest_1.5.tar.gz", repos=NULL, type="source")
install.packages("corpcor_1.6.10.tar.gz", repos=NULL, type="source")
install.packages("rARPACK_0.11-0.tar.gz", repos=NULL, type="source")
install.packages("ellipse_0.5.0.tar.gz", repos=NULL, type="source")
install.packages("mixOmics_6.28.0.tar.gz", repos=NULL, type="source")

library(mixOmics)


install.packages("ucminf_1.2.2.tar.gz", repos=NULL, type="source")
install.packages("ordinal_2023.12-4.1.tar.gz", repos=NULL, type="source")
install.packages("pan_1.9.tar.gz", repos=NULL, type="source")
install.packages("jomo_2.7-6.tar.gz", repos=NULL, type="source")
install.packages("mitml_0.4-5.tar.gz", repos=NULL, type="source")
install.packages("glmnet_4.1-8.tar.gz", repos=NULL, type="source")
install.packages("mice_3.16.0.tar.gz", repos=NULL, type="source")

library(mice)


install.packages("visNetwork_2.1.2.tar.gz", repos = NULL, type="source")
install.packages("modeltools_0.2-24.tar.gz", repos = NULL, type="source")
install.packages("flexmix_2.3-20.tar.gz", repos = NULL, type="source")

install.packages("prabclus_2.3-4.tar.gz", repos = NULL, type="source")
install.packages("diptest_0.77-1.tar.gz", repos = NULL, type="source")
install.packages("DEoptimR_1.1-3-1.tar.gz", repos = NULL, type="source")
install.packages("robustbase_0.99-4-1.tar.gz", repos = NULL, type="source")
install.packages("fpc_2.2-13.tar.gz", repos = NULL, type="source")
install.packages("aricode_1.0.3.tar.gz", repos = NULL, type="source")
install.packages("entropy_1.3.2.tar.gz", repos = NULL, type="source")
install.packages("terra_1.8-54.tar.gz", repos = NULL, type="source")
install.packages("raster_3.6-32.tar.gz", repos = NULL, type="source")

install.packages("proxy_0.4-27.tar.gz", repos = NULL, type="source")
install.packages("e1071_1.7-16.tar.gz", repos = NULL, type="source")
install.packages("classInt_0.4-11.tar.gz", repos = NULL, type="source")
install.packages("wk_0.9.4.tar.gz", repos = NULL, type="source")
install.packages("s2_1.1.9.tar.gz", repos = NULL, type="source")
install.packages("units_0.8-7.tar.gz", repos = NULL, type="source")

install.packages("sf_1.0-21.tar.gz", repos = NULL, type="source")
install.packages("sabre_0.4.3.tar.gz", repos = NULL, type="source")
install.packages("ggalluvial_0.12.5.tar.gz", repos = NULL, type="source")
install.packages("exactRankTests_0.8-35.tar.gz", repos = NULL, type="source")
install.packages("maxstat_0.7-26.tar.gz", repos = NULL, type="source")
install.packages("KMsurv_0.1-6.tar.gz", repos = NULL, type="source")
install.packages("km.ci_0.5-6.tar.gz", repos = NULL, type="source")
install.packages("survMisc_0.5.6.tar.gz", repos = NULL, type="source")
install.packages("litedown_0.5.tar.gz", repos = NULL, type="source")

install.packages("markdown_2.0.tar.gz", repos = NULL, type="source")
install.packages("gridtext_0.1.5.tar.gz", repos = NULL, type="source")
install.packages("ggtext_0.1.2.tar.gz", repos = NULL, type="source")

install.packages("survminer_0.5.0.tar.gz", repos = NULL, type="source")
install.packages("prettyGraphs_2.2.0.tar.gz", repos = NULL, type="source")
install.packages("ExPosition_2.11.0.tar.gz", repos = NULL, type="source")
install.packages("alluvial_0.1-2.tar.gz", repos = NULL, type="source")
install.packages("SNFtool_2.3.1.tar.gz", repos = NULL, type="source")
install.packages("shades_1.4.0.tar.gz", repos = NULL, type="source")
install.packages("ggfittext_0.10.2.tar.gz", repos = NULL, type="source")

library(SNFtool)
