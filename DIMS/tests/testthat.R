# Run all unit tests
library(testthat)
# enable snapshots
local_edition(3)

testthat::test_file("testthat/test_peak_finding_functions.R")
