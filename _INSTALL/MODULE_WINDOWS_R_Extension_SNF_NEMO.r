# #MANUAL INPUT
# args = as.list(c("BLymphocytes","SLE"))
# args[1] <-"mutated"
# args[2] <-"unmutated"
# args[3] <-"C:/Users/crist/Desktop/BiomiX2.5"
# 
# directory <- args[[3]]
# setwd(paste(directory,"_INSTALL",sep="/"))
# renv::load(paste(directory,"_INSTALL",sep="/"))

# library(masstools)
# library(metpath)
# detach("package:metpath", unload = TRUE)
# detach("package:masstools", unload = TRUE)


# Activate renv environment
cat("Activating renv environment...\n")
renv::activate()

# Show current library path
cat("📚 Installing into library:\n")
print(.libPaths())

setwd(paste(getwd(),"/", "Package_Windows/SNF_NEMO", sep=""))

chooseCRANmirror(48, ind = TRUE)

options(Ncpus = 6)



install.packages("visNetwork_2.1.2.zip", repos = NULL, type="source")
library(visNetwork)
install.packages("modeltools_0.2-24.zip", repos = NULL, type="source")
library(modeltools)
install.packages("flexmix_2.3-20.zip", repos = NULL, type="source")
library(flexmix)
install.packages("prabclus_2.3-4.zip", repos = NULL, type="source")
library(prabclus)
install.packages("diptest_0.77-1.zip", repos = NULL, type="source")
library(diptest)
install.packages("DEoptimR_1.1-3-1.zip", repos = NULL, type="source")
library(DEoptimR)
install.packages("robustbase_0.99-4-1.tar.gz", repos = NULL, type="source")
library(robustbase)
install.packages("kernlab_0.9-33.zip", repos = NULL, type="source")
library(kernlab)
install.packages("fpc_2.2-13.zip", repos = NULL, type="source")
library(fpc)
# 

install.packages("aricode_1.0.3.tar.gz", repos = NULL, type="source")
library(aricode)
install.packages("entropy_1.3.2.tar.gz", repos = NULL, type="source")
library(entropy)
install.packages("sp_1.6-1.tar.gz", repos = NULL, type="source")
library(sp)
# library(Rcpp)
# install.packages("Rcpp_1.0.14.tar.gz", repos = NULL, type="source")
# 
# 
# install.packages("terra_1.7-83.tar.gz", repos = NULL, type="source")
# 
# install.packages("terra_1.8-54.zip", repos = NULL, type="source")
# library(terra)
# install.packages("raster_3.6-32.zip", repos = NULL, type="source")
# library(raster)

install.packages("proxy_0.4-27.tar.gz", repos = NULL, type="source")
library(proxy)
install.packages("e1071_1.7-16.tar.gz", repos = NULL, type="source")
library(e1071)
install.packages("classInt_0.4-11.tar.gz", repos = NULL, type="source")
library(classInt)
install.packages("wk_0.9.4.tar.gz", repos = NULL, type="source")
library(wk)
install.packages("s2_1.1.9.zip", repos = NULL, type="source")
library(s2)
install.packages("units_0.8-7.tar.gz", repos = NULL, type="source")
library(units)
install.packages("sf_1.0-21.zip", repos = NULL, type="source")
library(sf)
# install.packages("sabre_0.4.3.zip", repos = NULL, type="source")
# install.packages("sabre_0.4.3.zip", repos = NULL, type="source")
# install.packages("sabre")
library(sabre)
install.packages("ggalluvial_0.12.5.zip", repos = NULL, type="source")
install.packages("exactRankTests_0.8-35.zip", repos = NULL, type="source")
install.packages("maxstat_0.7-26.zip", repos = NULL, type="source")
install.packages("KMsurv_0.1-6.zip", repos = NULL, type="source")
install.packages("km.ci_0.5-6.zip", repos = NULL, type="source")
install.packages("survMisc_0.5.6.zip", repos = NULL, type="source")
install.packages("litedown_0.3.tar.gz", repos = NULL, type="source")

install.packages("markdown_2.0.zip", repos = NULL, type="source")
install.packages("gridtext_0.1.5.zip", repos = NULL, type="source")
install.packages("ggtext_0.1.2.zip", repos = NULL, type="source")

install.packages("survminer_0.5.0.zip", repos = NULL, type="source")
install.packages("prettyGraphs_2.2.0.zip", repos = NULL, type="source")
install.packages("ExPosition_2.11.0.zip", repos = NULL, type="source")
install.packages("alluvial_0.1-2.zip", repos = NULL, type="source")
install.packages("SNFtool_2.3.1.zip", repos = NULL, type="source")
install.packages("shades_1.4.0.zip", repos = NULL, type="source")
install.packages("ggfittext_0.10.2.zip", repos = NULL, type="source")
# 
library(SNFtool)
install.packages("NEMO_0.1.0.tar.gz", repos = NULL, type="source")
library(NEMO)

print("MODULE SNF-NEMO INSTALLED")