# unit tests for PeakGrouping
# functions: find_peak_groups and annotate_peakgroups

# source all functions for PeakGrouping
source("../../preprocessing/peak_grouping_functions.R")

# test find_peak_groups
testthat::test_that("peak groups are correctly formed", {
  # create peak list to test on:
  samplenrs <- rep(c("C1", "P2", "P3"), 3)
  test_outlist_sorted_num <- matrix(0, nrow = 9, ncol = 5)
  colnames(test_outlist_sorted_num) <- c("mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")
  test_outlist_sorted_num[, 1] <- 300 + (1:9) / 1000
  test_outlist_sorted_num[, 5] <- 100 * (9:1)
  test_outlist_sorted <- as.data.frame(cbind(samplenr = samplenrs, test_outlist_sorted_num))
  test_outlist_sorted[, 2] <- as.numeric(test_outlist_sorted[, 2])
  test_outlist_sorted[, 6] <- as.numeric(test_outlist_sorted[, 6])
  test_mz_tolerance <- 0.0002
  test_sample_names <- unique(samplenrs)

  # test that first peak group is correctly formed
  expect_equal(as.numeric(find_peak_groups(test_outlist_sorted, test_mz_tolerance, test_sample_names)[1, c(5:7)]),
	       c(900, 0, 0), tolerance = 0.001, TRUE)
  # test column names
  expect_equal(colnames(find_peak_groups(test_outlist_sorted, test_mz_tolerance, test_sample_names)),
	       c("mzmed.pgrp", "mzmin.pgrp", "mzmax.pgrp", "nrsamples", unique(samplenrs)), TRUE)
  # test what happens if peak list is empty
  test_outlist_sorted_empty <- test_outlist_sorted[0, ]
  expect_equal(nrow(find_peak_groups(test_outlist_sorted_empty, test_mz_tolerance, test_sample_names)),
	       0, tolerance = 0.001, TRUE)
})

testthat::test_that("peak groups are correctly annotated", {
  # create peak group matrix
  test_peakgroup_matrix <- matrix(0, ncol = 7, nrow = 8)
  colnames(test_peakgroup_matrix) <- c("mzmed.pgrp", "mzmin.pgrp", "mzmax.pgrp", "nrsamples", "C1", "P2", "P3")
  test_peakgroup_matrix[, 1] <- 300 + (1:8) / 1000
  test_mz_tolerance <- 0.0001
  
  # create HMDB part
  test_hmdb_add_iso <- matrix(nrow = 4, ncol = 5)
  colnames(test_hmdb_add_iso) <- c("HMDB_ID_all", "sec_HMDB_ID", "CompoundName", "HMDB_name_all", "Mpos")
  rownames(test_hmdb_add_iso) <- c("HMDB1234561", "HMDB1234562_4", "iso_12345", "HMDB1234564")
  test_hmdb_add_iso[, 1] <- paste0("HMDB123456", 1:4)
  test_hmdb_add_iso[, 2] <- paste0("HMDB1234", 1:4)
  test_hmdb_add_iso[, 3] <- paste0("testcompound", 1:4)
  test_hmdb_add_iso[, 4] <- paste(test_hmdb_add_iso[, 3], test_hmdb_add_iso[, 3], sep = ";")
  test_hmdb_add_iso[, 5] <- 300 + (1:4) / 1000
  test_hmdb_add_iso_df <- as.data.frame(test_hmdb_add_iso)
  test_hmdb_add_iso_df[, 5] <- as.numeric(test_hmdb_add_iso_df[, 5])
  column_label <- "Mpos"
 
  # test that second peak group has been annotated
  expect_equal(as.character(annotate_peak_groups(test_peakgroup_matrix, test_hmdb_add_iso_df, column_label, test_mz_tolerance)[1, "assi_HMDB"]),
	       "testcompound1", TRUE)
  # test that third peak group has been annotated as isotope
  expect_equal(as.character(annotate_peak_groups(test_peakgroup_matrix, test_hmdb_add_iso_df, column_label, test_mz_tolerance)[3, "iso_HMDB"]),
	       "testcompound3", TRUE)
  # test what happens if hmdb part is empty
  test_hmdb_add_iso_empty <- test_hmdb_add_iso_df[0, ]
  expect_equal(as.logical(annotate_peak_groups(test_peakgroup_matrix, test_hmdb_add_iso_empty, column_label, test_mz_tolerance)[, "assi_HMDB"]),
	       rep(NA, 8), TRUE)
})
