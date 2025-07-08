# unit tests for PeakFinding functions:

source("../preprocessing/peak_finding_functions.R")

# test search_regions_of_interest function:
test_that("regions of interest are correctly found for a single peak", {
  # create intensities to test on
  test_stdev <- 0.00006
  test_mass_vector <- rep(70, 10) + 1:10 * test_stdev
  test_ints_range <- dnorm(test_mass_vector, mean = 70.0003, sd = test_stdev)
  test_ints_fullrange <- as.data.frame(cbind(mz = test_mass_vector, 
                                             int = test_ints_range))
  
  # expected output
  expected_output <- matrix(c(1, 10, 10), nrow = 1, ncol = 3L)
  colnames(expected_output) <- c("from", "to", "length")
  rownames(expected_output) <- "to"

  # tests
  expect_type(search_regions_of_interest(test_ints_fullrange), "double")
  expect_equal(nrow(search_regions_of_interest(test_ints_fullrange)), 1)
  expect_equal(search_regions_of_interest(test_ints_fullrange), expected_output, tolerance = 0.000001, TRUE)
})
test_that("regions of interest are correctly found for two peaks", {
  # create intensities to test on
  test_stdev <- 0.00006
  test_mass_vector <- rep(70, 20) + 1:20 * test_stdev
  test_ints_range <- dnorm(test_mass_vector, mean = 70.0003, sd = test_stdev) + 
                     dnorm(test_mass_vector, mean = 70.0005, sd = test_stdev)
  test_ints_fullrange <- as.data.frame(cbind(mz = test_mass_vector, 
                                             int = test_ints_range))
  
  # expected output
  expected_output <- as.data.frame(matrix(c(1, 8, 8, 20, 8, 13), nrow = 2, ncol = 3))
  colnames(expected_output) <- c("from", "to", "length")
  rownames(expected_output) <- as.character(c(1, 2))
  
  # tests
  expect_type(search_regions_of_interest(test_ints_fullrange), "list")
  expect_equal(nrow(search_regions_of_interest(test_ints_fullrange)), 2)
  expect_equal(search_regions_of_interest(test_ints_fullrange)[, "length"], as.numeric(expected_output[, "length"]))
  expect_equal(search_regions_of_interest(test_ints_fullrange), expected_output, check.attributes = FALSE)
})

# test integrate_peaks function:
test_that("peaks are correctly integrated", {
  # create peak info to test on:
  test_stdev <- 0.00006
  test_mass_vector <- rep(70, 20) + 1:20 * test_stdev
  test_ints_range <- dnorm(test_mass_vector, mean = 70.0003, sd = test_stdev) + 
                     dnorm(test_mass_vector, mean = 70.0005, sd = test_stdev)
  test_ints_fullrange <- as.data.frame(cbind(mz = test_mass_vector, 
                                             int = test_ints_range))
  
  # create regions of interest to test on:
  test_regions_of_interest <- as.data.frame(matrix(c(1, 8, 8, 20, 8, 13), nrow = 2, ncol = 3L))
  colnames(test_regions_of_interest) <- c("from", "to", "length")
  rownames(test_regions_of_interest) <- as.character(c(1, 2))
  
  # set other parameters
  test_resol <- 140000
  test_peak_thresh <- 2000
  
  # expected output
  expected_output <- as.data.frame(matrix(c(70.000357, 70.000521, 0.5148094, 0.3016569, 70.00006, 70.00048, 
                                            70.00048, 70.00120, 15526.826791, 11950.90570), nrow = 2, ncol = 5))
  colnames(expected_output) <- c("mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")
  
  expect_type(integrate_peaks(test_ints_fullrange, test_regions_of_interest, test_resol, test_peak_thresh), "double")
  expect_equal(nrow(integrate_peaks(test_ints_fullrange, test_regions_of_interest, test_resol, test_peak_thresh)), 2)
  expect_equal(integrate_peaks(test_ints_fullrange, test_regions_of_interest, test_resol, test_peak_thresh)[, "mzmed.pkt"], 
               expected_output[, "mzmed.pkt"], tolerance = 0.0001)
  expect_equal(integrate_peaks(test_ints_fullrange, test_regions_of_interest, test_resol, test_peak_thresh)[, "height.pkt"], 
               expected_output[, "height.pkt"], tolerance = 0.0001)
})

# test get_fwhm function
test_that("fwhm is correctly calculated", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_mz_max <- max(test_mass_vector)
  test_resol <- 140000

  expect_type(get_fwhm(test_mz_max, test_resol), "double")
  expect_equal(get_fwhm(test_mz_max, test_resol), 0.000295865, tolerance = 0.000001, TRUE)
})

