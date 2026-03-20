# unit test for AveragePeaks

# source all functions for PeakGrouping
source("../../preprocessing/average_peaks_functions.R")

# test find_peak_groups
testthat::test_that("peaks are correctly averaged", {
  # create peak list to test on:
  samplenrs <- rep(c("repl_001", "repl_002", "repl_003"), 3)
  test_peaklist_sorted_num <- matrix(0, nrow = 9, ncol = 5)
  colnames(test_peaklist_sorted_num) <- c("mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")
  test_peaklist_sorted_num[, 1] <- 300 + (1:9) / 10000
  test_peaklist_sorted_num[, 5] <- 100 * (9:1)
  test_peaklist_sorted <- as.data.frame(cbind(samplenr = samplenrs, test_peaklist_sorted_num))
  test_peaklist_sorted[, 2] <- as.numeric(test_peaklist_sorted[, 2])
  test_peaklist_sorted[, 6] <- as.numeric(test_peaklist_sorted[, 6])
  test_sample_name <- "P001"
  test_empty_peaklist <- test_peaklist_sorted[0, ]

  # test that first peak is correctly averaged
  expect_equal(as.numeric(average_peaks_per_sample(test_peaklist_sorted, test_sample_name)[1, 6]), 600, tolerance = 0.001, TRUE)
  # test number of rows
  expect_equal(nrow(average_peaks_per_sample(test_peaklist_sorted, test_sample_name)), 2)
  # test column names
  expect_equal(colnames(average_peaks_per_sample(test_peaklist_sorted, test_sample_name)),
               c("samplenr", "mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt"), TRUE)
  # test what happens when peak list is empty
  expect_error(average_peaks_per_sample(test_empty_peaklist, test_sample_name))
})

