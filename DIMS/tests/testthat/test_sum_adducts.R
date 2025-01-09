# test function sum_adducts

# create peakgroup_list to test on; diagnostics setting
test_peakgroup_list <- data.frame(matrix(NA, nrow = 4, ncol = 23))
colnames(test_peakgroup_list) <- c("mzmed.pgrp", "nrsamples", "ppmdev", "assi_HMDB", "all_hmdb_names",
                                   "iso_HMDB", "HMDB_code", "all_hmdb_ids", "sec_hmdb_ids", "theormz_HMDB",
                                   "C101.1", "C102.1", "P2.1", "P3.1",
                                   "avg.int", "assi_noise", "theormz_noise", "avg.ctrls", "sd.ctrls",
                                   "C101.1_Zscore", "C102.1_Zscore", "P2.1_Zscore", "P3.1_Zscore")
# fill relevant columns
test_peakgroup_list[ , c(1)] <- 300 + runif(4)
test_peakgroup_list[ , c(2, 3)] <- runif(8)
test_peakgroup_list[ , "HMDB_code"] <- c("HMDB1234567", "HMDB1234567_1", "HMDB1234567_2", "HMDB1234567_7")
test_peakgroup_list[ , "all_hmdb_names"] <- paste(test_peakgroup_list[ , "assi_HMDB"], 
                                                  test_peakgroup_list[ , "assi_HMDB"], sep=";")
test_peakgroup_list[ , grep("C", colnames(test_peakgroup_list))] <- 1000*(1:16)
test_peakgroup_list[ , grep("P", colnames(test_peakgroup_list))] <- 10000*(1:16)

# create HMDB part object 
test_hmdb_main_part <- matrix(NA, nrow = 2, ncol = 8)
colnames(test_hmdb_main_part) <- c("HMDB_ID_all", "sec_HMDB_ID", "CompoundName", "HMDB_name_all",
                                   "Composition", "MNeutral", "MNeg", "Mpos")
rownames(test_hmdb_main_part) <- c("HMDB1234567", "HMDB7654321")
test_hmdb_main_part[, 1] <- c("HMDB1234567;HMDB0000567", "HMDB7654321;HMDB0000321")

test_that("input is a directory", expect_equal(dir.exists("./testthat"), TRUE))
# use the function and test the output
test_sum_adduct_output <- sum_adducts(test_peakgroup_list, test_hmdb_main_part, c(1, 2), 1)
test_that("adduct sums are correctly calculated", expect_equal(as.vector(test_sum_adduct_output[1, c(1:4)]), c("6000", "18000", "60000", "180000"), TRUE))

