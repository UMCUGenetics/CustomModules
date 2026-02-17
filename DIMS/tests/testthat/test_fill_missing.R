# unit tests for FillMissing
# function: fill_missing_intensities
source("../../preprocessing/fill_missing_functions.R")

# test fill_missing_intensities
testthat::test_that("missing values are corretly filled with random values", {
  # It's necessary to copy/symlink the files to the current location for the combine_sum_adducts_parts function
  test_files <- list.files("fixtures/", "test_peakgroup_list", full.names = TRUE)
  file.symlink(file.path(test_files), getwd())
  
  # create replication pattern of technical replicates
  test_repl_pattern <- c(list(1), list(2), list(3), list(4))
  names(test_repl_pattern) <- c("C101.1", "C102.1", "P2.1", "P3.1")
  
  # read in peakgroup_list, set intensities for patient columns to zero
  test_peakgroup_list <- read.table("./test_peakgroup_list.txt", sep= "\t")
  test_peakgroup_list[, grep("P", colnames(test_peakgroup_list))] <- 0
  
  test_thresh <- 2000

  # create a large peak group list to test for negative values
  test_large_peakgroup_list <- rbind(test_peakgroup_list, test_peakgroup_list)
  for (index in 1:15) {
    test_large_peakgroup_list <- rbind(test_large_peakgroup_list, test_large_peakgroup_list)
  }
  # for the sake of time, leave only one intensity column with zeros
  test_large_peakgroup_list$P2.1 <- 1

  expect_equal(round(fill_missing_intensities(test_peakgroup_list, test_repl_pattern, test_thresh, not_random = TRUE)$P2.1),
               c(1944, 1977, 2156, 2007), TRUE, tolerance = 0.1)
  # fill_missing_intensities should not produce any negative values, even if a large quantity of numbers are filled in
  expect_gt(min(fill_missing_intensities(test_large_peakgroup_list, test_repl_pattern, test_thresh, not_random = FALSE)$P3.1),
            0, TRUE)

})

