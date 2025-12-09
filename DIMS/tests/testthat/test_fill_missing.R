# unit tests for FillMissing
# function: fill_missing_intensities
source("../../preprocessing/fill_missing_functions.R")

# test fill_missing_intensities
testthat::test_that("missing values are corretly filled with random values", {
  # create peakgroup_list to test on in diagnostics setting
  test_peakgroup_list <- data.frame(matrix(NA, nrow = 4, ncol = 23))
  colnames(test_peakgroup_list) <- c("mzmed.pgrp", "nrsamples", "ppmdev", "assi_HMDB", "all_hmdb_names",
                                     "iso_HMDB", "HMDB_code", "all_hmdb_ids", "sec_hmdb_ids", "theormz_HMDB",
                                     "C101.1", "C102.1", "P2.1", "P3.1",
                                     "avg.int", "assi_noise", "theormz_noise", "avg.ctrls", "sd.ctrls",
                                     "C101.1_Zscore", "C102.1_Zscore", "P2.1_Zscore", "P3.1_Zscore")
  test_peakgroup_list[, c(1)] <- 300 + runif(4)
  test_peakgroup_list[, c(2, 3)] <- runif(8)
  test_peakgroup_list[, "HMDB_code"] <- c("HMDB1234567", "HMDB1234567_1", "HMDB1234567_2", "HMDB1234567_7")
  test_peakgroup_list[, "all_hmdb_ids"] <- paste(test_peakgroup_list[, "HMDB_code"],
                                                 test_peakgroup_list[, "HMDB_code"], sep = ";")
  test_peakgroup_list[, "all_hmdb_names"] <- paste(test_peakgroup_list[, "assi_HMDB"],
                                                   test_peakgroup_list[, "assi_HMDB"], sep = ";")
  test_peakgroup_list[, grep("C", colnames(test_peakgroup_list))] <- 1000 * (1:16)
  test_peakgroup_list[, grep("P", colnames(test_peakgroup_list))] <- 0
  test_repl_pattern <- c(list(1), list(2), list(3), list(4))
  names(test_repl_pattern) <- c("C101.1", "C102.1", "P2.1", "P3.1")
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
  start.time <- Sys.time()
  expect_gt(min(fill_missing_intensities(test_large_peakgroup_list, test_repl_pattern, test_thresh, not_random = FALSE)$P3.1),
            0, TRUE)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken

})
