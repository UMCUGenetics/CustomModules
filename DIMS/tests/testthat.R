# Run all unit tests on sum_adducts
library(testthat)

# source all functions for running Outrider
source("../Utils/sum_adducts.R")

testthat::test_file("testthat/test_sum_intensities_adducts.R")

