# Run all unit tests
library(testthat)
# enable snapshots
local_edition(3)

testthat::test_file("testthat/test_check_qc.R")
testthat::test_file("testthat/test_collect_sum_adducts.R")
testthat::test_file("testthat/test_generate_excel.R")
