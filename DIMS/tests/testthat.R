# Run all unit tests
library(testthat)
library(withr)
# enable snapshots
local_edition(3)

Sys.setenv(NOT_CRAN = "true")
testthat::test_dir("tests/testthat")
