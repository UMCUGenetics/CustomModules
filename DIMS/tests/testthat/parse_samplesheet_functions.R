# unit tests for ParseSamplesheet
# function: generate_repl_pattern

# source all functions for ParseSamplesheet
source("../../preprocessing/parse_samplesheet_functions.R")

# test generate_repl_pattern
testthat::test_that("replication pattern is correctly generated", {
  # create sample sheet tot test on:
  test_file_names <- paste0(rep("RES_20260101_", 6), sprintf("%03d", 1:6))
  test_sample_names <- sort(rep(c("C1", "P2", "P3"), 2))
  test_sample_sheet <- as.data.frame(cbind(File_Name = test_file_names, Sample_Name = test_sample_names))

  # test that a list of length 3 is generated
  expect_length(generate_repl_pattern(test_sample_sheet), 3)
  # test list names
  expect_equal(names(generate_repl_pattern(test_sample_sheet)), unique(test_sample_names), TRUE)
  
  # test what happens if any sample name is used twice
  test_sample_names <- gsub("P3", "P2", test_sample_names)
  test_sample_sheet <- as.data.frame(cbind(File_Name = test_file_names, Sample_Name = test_sample_names))
  expect_length(generate_repl_pattern(test_sample_sheet), 2)
  expect_length(generate_repl_pattern(test_sample_sheet)$P2, 4)
})
