# unit tests for CollectFilled
# functions: collapse_information, merge_duplicate_rows, calculate_zscores_peakgrouplist,
#            calculate_ppm_deviation, order_columns_peakgrouplist
source("../../preprocessing/collect_filled_functions.R")

testthat::test_that("Values for duplicate rows are correctly collapsed", {
  test_matrix <- matrix(letters[1:8], nrow = 2, ncol = 4)
  colnames(test_matrix) <- paste0("column", 1:4)

  expect_equal(collapse_information("column1", test_matrix, c(1,2)), "a;b", TRUE)
})

testthat::test_that("Duplicate rows in a peak group list are correctly merged", {
  # Copy/symlink the files to the current location for the function
  test_files <- list.files("fixtures/", "test_peakgroup_list", full.names = TRUE)
  file.symlink(file.path(test_files), getwd())

  # read in peakgroup_list, create duplicate row 
  test_peakgroup_list <- read.table("./test_peakgroup_list.txt", sep = "\t")
  test_peakgroup_list_dup <- test_peakgroup_list[c(1, 2, 2, 3), ]

  # after merging duplicate rows, the test peak group list should have 3 rows
  expect_equal(nrow(merge_duplicate_rows(test_peakgroup_list_dup)), 3, TRUE, tolerance = 0.001)
  expect_equal(merge_duplicate_rows(test_peakgroup_list_dup)[3, "all_hmdb_ids"], 
               paste(test_peakgroup_list_dup[2, "all_hmdb_ids"], test_peakgroup_list_dup[3, "all_hmdb_ids"], sep = ";"), 
               TRUE)
})

testthat::test_that("Z-scores are correctly calculated in CollectFilled", {
  # read in peakgroup_list; remove avg.int, std.int, noise and Zscore columns
  test_peakgroup_list <- read.table("./test_peakgroup_list.txt", sep = "\t")
  # remove avg.ctrls, sd.ctrls and Z-score columns
  test_peakgroup_list_noz <- test_peakgroup_list[ , -grep("avg.ctrls", colnames(test_peakgroup_list))]
  test_peakgroup_list_noz <- test_peakgroup_list_noz[ , -grep("sd.ctrls", colnames(test_peakgroup_list_noz))]
  test_peakgroup_list_noz <- test_peakgroup_list_noz[ , -grep("_Zscore", colnames(test_peakgroup_list_noz))]
 
  # after calculate_zscores_peakgrouplist, there should be 4 columns with _Zscore in the name
  expect_equal(length(grep("_Zscore", colnames(calculate_zscores_peakgrouplist(test_peakgroup_list_noz)))), 4, TRUE, tolerance = 0.001)

  # after calculate_zscores_peakgrouplist, the 4 columns with _Zscore in the name should be filled non-zero
  expect_equal(calculate_zscores_peakgrouplist(test_peakgroup_list_noz)$C101.1_Zscore[1], -0.7071, TRUE, tolerance = 0.00001)
  expect_equal(calculate_zscores_peakgrouplist(test_peakgroup_list_noz)$P2.1_Zscore[4], 12.0208, TRUE, tolerance = 0.00001)
 
})

testthat::test_that("ppm deviation values are correctly calculated in CollectFilled", {
  # read in peakgroup_list
  test_peakgroup_list <- read.table("./test_peakgroup_list.txt", sep = "\t")

  # store ppm deviation values
  test_ppm_values <- test_peakgroup_list$ppmdev

  # after calculate_ppm_deviation, there ppm values in the new column should approximate the old ones
  expect_equal(calculate_ppm_deviation(test_peakgroup_list)$ppmdev, test_ppm_values, TRUE, tolerance = 0.001)

  # calculate_ppm_deviation should give NA if there is no value for theormz_HMDB
  test_peakgroup_list$theormz_HMDB[1] <- NA
  expect_identical(is.na(calculate_ppm_deviation(test_peakgroup_list)$ppmdev[1]), TRUE)
 
})

testthat::test_that("columns in peak group list are corretly sorted", {
  # read in peakgroup_list
  test_peakgroup_list <- read.table("./test_peakgroup_list.txt", sep = "\t")
  # original order of columns
  original_column_order <- colnames(test_peakgroup_list)
  # after ordering, column names should be re-ordered
  test_column_order <- original_column_order[c(1, 2, 7:14, 21, 3:6, 15:20)]
 
  expect_identical(colnames(order_columns_peakgrouplist(test_peakgroup_list)), test_column_order)

  # Remove copied/symlinked files
  files_remove <- list.files("./", "test_peakgroup_list.txt", full.names = TRUE)
  file.remove(files_remove)
 
})

