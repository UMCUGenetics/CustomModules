# unit tests for GenerateBreaks
# functions: get_trim_parameters and get_breaks_for_bins

# source all functions for GenerateBreaks
source("../../preprocessing/generate_breaks_functions.R")

# test get_trim_parameters
testthat::test_that("trim parameters are correctly calculated", {
  # create list of scan times to test on:
  test_scantimes <- seq(0, 100, length = 168)
  test_polarities <- c(rep("positive", 84), rep("negative", 84))
  test_trim <- 0.1

  # test that the function produces no output except trim_params.RData file
  expect_silent(get_trim_parameters(test_scantimes, test_polarities, test_trim))
  expect_true(file.exists("./trim_params.RData"))

  # test the parameters in the output RData file
  load("./trim_params.RData")
  test_trim_left_pos  <- trim_left_pos
  test_trim_left_neg  <- trim_left_neg
  test_trim_right_pos <- trim_right_pos
  test_trim_right_neg <- trim_right_neg
  rm(trim_left_pos, trim_left_neg, trim_right_pos, trim_right_neg)

  # load previously stored values from fixtures
  load("fixtures/test_trim_params.RData")
  expect_equal(test_trim_left_pos,  trim_left_pos)
  expect_equal(test_trim_left_neg,  trim_left_neg)
  expect_equal(test_trim_right_pos, trim_right_pos)
  expect_equal(test_trim_right_neg, trim_right_neg)
})

# test get_breaks_for_bins
testthat::test_that("breaks are correctly calculated", {
  # create list of scan times to test on:
  test_mzrange <- c(300, 400)
  test_resol <- 140000
  
  # test that the function produces no output except breaks.fwhm.RData file
  expect_silent(get_breaks_for_bins(test_mzrange, test_resol))
  expect_true(file.exists("./breaks.fwhm.RData"))
  
  # test the vectors in the output RData file
  load("./breaks.fwhm.RData")
  test_breaks_fwhm <- breaks_fwhm
  test_breaks_fwhm_avg <- breaks_fwhm_avg
  rm(breaks_fwhm, breaks_fwhm_avg)

  # load breaks from fixtures and compare vectors
  load("fixtures/breaks.fwhm.RData")
  expect_equal(test_breaks_fwhm, breaks_fwhm)
  expect_equal(test_breaks_fwhm_avg, breaks_fwhm_avg)
})
