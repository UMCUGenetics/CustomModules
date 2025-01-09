# Run all unit tests on sum_adducts
library(testthat)

# source all functions for running Outrider
source("/Users/mraves2/Development/DIMS_refactor_SumAdducts/DIMS/CustomModules/DIMS/Utils/sum_adducts.R")

# temporary (files are local):
# setwd("/Users/mraves2/Genomics/RNAseq_outrider/tests")

# set working directory to location of testthat.R file
location_current_file <- rstudioapi::getSourceEditorContext()$path
path_to_current_file <- strsplit(location_current_file, "testthat.R")[[1]][1]
setwd(path_to_current_file)

testthat::test_file("testthat/test_sum_adducts.R")

