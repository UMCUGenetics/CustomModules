# unit tests for CollectSumAdducts
# functions: combine_sum_adduct_parts, combine_scanmodes_intensities
source("../../preprocessing/collect_sum_adducts_functions.R")


testthat::test_that("Combining summed adduct parts", {
  # It's necessary to copy/symlink the files to the current location for the combine_sum_adducts_parts function
  test_files <- list.files("test_data/", "SummedAdducts_test.RData", full.names = TRUE)
  file.symlink(file.path(test_files), getwd())

  expect_equal(rownames(combine_sum_adduct_parts("positive")),
               c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB005", "HMDB006", "HMDB007", "HMDB008",
                 "HMDB009", "HMDB010", "HMDB012", "HMDB014"))
  expect_equal(colnames(combine_sum_adduct_parts("positive")),
               c("C101.1", "C102.1", "P2.1", "P3.1", "HMDB_name", "HMDB_ID_all", "sec_HMDB_ID", "HMDB_name_all"))

  # Remove copied/symlinked files
  files_remove <- list.files("./", "SummedAdducts_test.RData", full.names = TRUE)
  file.remove(files_remove)
})


testthat::test_that("Combining intensities of positive and negative scanmodes", {
  test_outlist_pos <- data.frame(
    C101.1 = c(100, 200, 300, 400),
    C102.1 = c(125, 225, 325, 425),
    P2.1 = c(150, 250, 350, 450),
    P3.1 = c(175, 275, 375, 475),
    HMDB_name = c("metab_1", "metab_2", "metab_3", "metab_4"),
    HMDB_ID_all = c("HMDB001;HMDB011", "HMDB002", "HMDB003;HMDB013", "HMDB004"),
    sec_HMDB_ID = c("HMDB1;HMDB11", "", "HMDB3;HMDB13", "HMDB4"),
    HMDB_name_all = c("metab_1;metab_11", "metab_2", "metab_3;metab_13", "metab_4")
  )
  rownames(test_outlist_pos) <- c("HMDB001", "HMDB002", "HMDB003", "HMDB004")
  test_outlist_neg <- data.frame(
    C101.1 = c(100, 200, 300, 400),
    C102.1 = c(125, 225, 325, 425),
    P2.1 = c(150, 250, 350, 450),
    P4.1 = c(175, 275, 375, 475),
    HMDB_name = c("metab_1", "metab_2", "metab_3", "metab_5"),
    HMDB_ID_all = c("HMDB001;HMDB011", "HMDB002", "HMDB003;HMDB013", "HMDB005"),
    sec_HMDB_ID = c("HMDB1;HMDB11", "", "HMDB3;HMDB13", "HMDB5"),
    HMDB_name_all = c("metab_1;metab_11", "metab_2", "metab_3;metab_13", "metab_5")
  )
  rownames(test_outlist_neg) <- c("HMDB001", "HMDB002", "HMDB003", "HMDB005")

  # Check that P3.1 and P4.1 are not present
  expect_equal(colnames(combine_scanmodes_intensities(test_outlist_pos, test_outlist_neg)),
               c("C101.1", "C102.1", "P2.1", "HMDB_name", "HMDB_name_all", "HMDB_ID_all", "sec_HMDB_ID"))
  # Check if combined intensities are correct
  expect_equal(combine_scanmodes_intensities(test_outlist_pos, test_outlist_neg)$C101.1, c(200, 400, 600, 400, 400))
  # Check if all metabolites are present in the rownames
  expect_equal(rownames(combine_scanmodes_intensities(test_outlist_pos, test_outlist_neg)),
               c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB005"))
})
