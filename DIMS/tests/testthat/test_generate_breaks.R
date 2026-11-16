# Unit tests for GenerateBreaks
# functions: get_trim_parameters and get_breaks_for_bins

# Source all functions for GenerateBreaks
source("../../preprocessing/generate_breaks_functions.R")

# Test get_trim_parameters
testthat::test_that("trim parameters are correctly calculated", {
  # Create list of scan times to test on:
  test_scantimes <- seq(0, 100, length = 168)
  test_polarities <- c(rep("positive", 84), rep("negative", 84))
  test_trim <- 0.1

  # Test that the function produces no output except trim_params.RData file
  expect_silent(get_trim_parameters(test_scantimes, test_polarities, test_trim))
  expect_true(file.exists("./trim_params.RData"))

  # Test the parameters in the output RData file
  load("./trim_params.RData")
  test_trim_left_pos  <- trim_left_pos
  test_trim_left_neg  <- trim_left_neg
  test_trim_right_pos <- trim_right_pos
  test_trim_right_neg <- trim_right_neg
  rm(trim_left_pos, trim_left_neg, trim_right_pos, trim_right_neg)

  # Load previously stored values from fixtures
  load("fixtures/test_trim_params.RData")
  expect_equal(test_trim_left_pos,  trim_left_pos)
  expect_equal(test_trim_left_neg,  trim_left_neg)
  expect_equal(test_trim_right_pos, trim_right_pos)
  expect_equal(test_trim_right_neg, trim_right_neg)
})

# Test get_breaks_for_bins
testthat::test_that("breaks are correctly calculated", {
  # Create list of scan times to test on:
  test_mzrange <- c(300, 400)
  test_resol <- 140000
  
  # Test that the function produces no output except breaks.fwhm.RData file
  expect_silent(get_breaks_for_bins(test_mzrange, test_resol))
  expect_true(file.exists("./breaks.fwhm.RData"))
  
  # Test the vectors in the output RData file
  load("./breaks.fwhm.RData")
  test_breaks_fwhm <- breaks_fwhm
  test_breaks_fwhm_avg <- breaks_fwhm_avg
  rm(breaks_fwhm, breaks_fwhm_avg)

  # Load breaks from fixtures and compare vectors
  load("fixtures/test_breaks.fwhm.RData")
  expect_equal(test_breaks_fwhm, breaks_fwhm)
  expect_equal(test_breaks_fwhm_avg, breaks_fwhm_avg)
})
