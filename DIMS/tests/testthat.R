# Run all unit tests
library(testthat)
library(withr)
# enable snapshots
local_edition(3)

Sys.setenv(NOT_CRAN = "true")
results_tests <- test_dir("tests/testthat", reporter = "summary")

testthat::snapshot_review('check_qc/')

if (any(vapply(results_tests, function(x) any(x$failed > 0), logical(1)))) {
  quit(status = 1)
}
