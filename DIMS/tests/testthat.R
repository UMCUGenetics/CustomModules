# Run all unit tests on sum_adducts
library(testthat)

# source all functions for running Outrider
source("../Utils/sum_adducts.R")
source("../preprocessing/peak_grouping_functions.R")

testthat::test_file("testthat/test_sum_adducts.R")
testthat::test_file("testthat/test_find_peak_grouping.R")

