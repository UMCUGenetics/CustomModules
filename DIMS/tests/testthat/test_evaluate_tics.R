# unit test for EvaluateTics

# source all functions for PeakGrouping
source("../../preprocessing/evaluate_tics_functions.R")

# test find_bad_replicates
testthat::test_that("TICS are correctly accepted or rejected", {
  # It's necessary to copy/symlink the files to the current location for the combine_sum_adducts_parts function
  # local: setwd("~/Development/DIMS_refactor_PeakFinding_codereview/CustomModules/DIMS/tests/testthat")
  test_files <- list.files("fixtures/", "test_evaluate_tics", full.names = TRUE)
  file.symlink(file.path(test_files), getwd())
  
  # create replication pattern to test on:
  technical_replicates <- paste0("test_evaluate_tics_file", 1:3)
  test_repl_pattern <- list(technical_replicates)
  names(test_repl_pattern) <- "sample1"
  test_thresh2remove <- 10^9
  
  # test that output has two entries
  expect_equal(length(find_bad_replicates(test_repl_pattern, test_thresh2remove)), 2)
  # test that first technical replicate is removed in positive scan mode
  expect_equal(find_bad_replicates(test_repl_pattern, test_thresh2remove)$pos, "test_evaluate_tics_file1", TRUE)
  # test that third technical replicate is removed in negative scan mode
  expect_equal(find_bad_replicates(test_repl_pattern, test_thresh2remove)$neg, "test_evaluate_tics_file3", TRUE)
  # test that output files are generated
  expect_equal(sum(grepl("miss_infusions_", list.files("./"))), 2)
  
  # Remove symlinked files
  files_remove <- list.files("./", "SummedAdducts_test.RData", full.names = TRUE)
  file.remove(files_remove)
  
})

# test remove_from_repl_pattern
testthat::test_that("technical replicates are correctly removed from replication pattern", {
  # create replication pattern to test on:
  technical_replicates <- paste0("test_evaluate_tics_file", 1:3)
  test_repl_pattern <- list(technical_replicates)
  names(test_repl_pattern) <- "sample1"
  test_bad_samples <- "test_evaluate_tics_file2"
  test_nr_replicates <- 3
  
  # test that the output contains 1 sample
  expect_equal(length(remove_from_repl_pattern(test_bad_samples, test_repl_pattern, test_nr_replicates)), 1)
  # test that the output for the sample contains 2 technical replicates
  expect_equal(length(remove_from_repl_pattern(test_bad_samples, test_repl_pattern, test_nr_replicates)$sample1), 2)
  # test that the second technical replicate has been removed
  expect_false(unique(grepl(test_bad_samples, remove_from_repl_pattern(test_bad_samples, test_repl_pattern, test_nr_replicates)$sample1)))
})

# test get_overview_tech_reps
testthat::test_that("overview of technical replicates is correctly created", {
  # create replication pattern to test on:
  technical_replicates <- paste0("test_evaluate_tics_file", 1:3)
  test_repl_pattern_filtered <- list(technical_replicates)
  names(test_repl_pattern_filtered) <- "sample1"
  test_scanmode <- "positive"

  # test that overview is correctly created
  expect_equal(get_overview_tech_reps(test_repl_pattern_filtered, test_scanmode)[ ,1], "sample1")
  expect_equal(get_overview_tech_reps(test_repl_pattern_filtered, test_scanmode)[ ,3], "positive")
  expect_true(get_overview_tech_reps(test_repl_pattern_filtered, test_scanmode)[ ,2] == 
               paste0(technical_replicates, collapse = ";"), TRUE)
  
})

