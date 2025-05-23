# Run all unit tests
library(testthat)
library(withr)
# enable snapshots
local_edition(3)

# testthat::test_file("tests/testthat/test_check_qc.R")
# testthat::test_file("tests/testthat/test_collect_sum_adducts.R")
# testthat::test_file("tests/testthat/test_generate_excel.R")
# testthat::test_file("tests/testthat/test_sum_intensities_adducts.R")

results_tests <- test_dir("tests/testthat", reporter = "summary")

if (any(vapply(results_tests, function(x) anyx$failed > 0, logical(1)))) {
  quit(status = 1)
}
