# unit tests for PeakFinding functions: 

source("../preprocessing/peak_finding_functions.R")

# test fit_quality function:
test_that("fit quality is correctly calculated", {
  # initialize
  test_resol <- 140000
  test_scale_factor <- 2
  test_sigma <- 0.0001
  test_sum_fit <- NULL

  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_int_vector <- c(9000, 16000, 19000, 15000, 6000)
  test_mz_max <- max(test_mass_vector)

  expect_type(get_fit_quality(test_mass_vector, test_int_vector, test_mz_max, test_resol, test_scale_factor, test_sigma, test_sum_fit), "list")
  expect_equal(get_fit_quality(test_mass_vector, test_int_vector, test_mz_max, test_resol, test_scale_factor, test_sigma, test_sum_fit)$fq_new, 2.426, tolerance = 0.001, TRUE)
})

# test get_stdev:
test_that("standard deviation is correctly calculated", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_int_vector <- c(9000, 16000, 19000, 15000, 6000)
  test_resol <- 140000

  expect_type(get_stdev(test_mass_vector, test_int_vector, test_resol), "double")
  expect_equal(get_stdev(test_mass_vector, test_int_vector, test_resol), 0.000126, tolerance = 0.00001, TRUE)
})

# test estimate_area
test_that("area under the curve is correctly calculated", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_mz_max <- max(test_mass_vector)
  test_resol <- 140000
  test_scale_factor <- 2
  test_sigma <- 0.0001

  expect_type(estimate_area(test_mz_max, test_resol, test_scale_factor, test_sigma), "double")
  expect_equal(estimate_area(test_mz_max, test_resol, test_scale_factor, test_sigma), 1673.061, tolerance = 0.001, TRUE)
})

# test get_fwhm
test_that("fwhm is correctly calculated", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_mz_max <- max(test_mass_vector)
  test_resol <- 140000

  expect_type(get_fwhm(test_mz_max, test_resol), "double")
  expect_equal(get_fwhm(test_mz_max, test_resol), 0.000295865, tolerance = 0.000001, TRUE)
})

# test get_fwhm for NA
test_that("fwhm for mass NA gives standard fwhm", {
  test_resol <- 140000

  expect_type(get_fwhm(NA, test_resol), "double")
  expect_equal(get_fwhm(NA, test_resol), 0.001428571, tolerance = 0.000001, TRUE)
})

# test within_ppm
test_that("two peaks are within 5 ppm", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_mz_max <- max(test_mass_vector)
  test_scale_factor <- 2
  test_sigma <- 0.0001
  test_area <- 20000
  test_mass_vector_eq <- seq(min(test_mass_vector), max(test_mass_vector), length = 100)
  test_ppm <- 5
  test_resol <- 140000

  expect_type(within_ppm(test_mz_max, test_scale_factor, test_sigma, test_area, test_mass_vector_eq, test_mass_vector, test_ppm, test_resol), "list")
  expect_equal(within_ppm(test_mz_max, test_scale_factor, test_sigma, test_area, test_mass_vector_eq, test_mass_vector, test_ppm, test_resol)$mean, 70.00962, tolerance = 0.0001, TRUE)
})

# test sum_curves
test_that("two curves are correctly summed", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_mass_vector_eq <- seq(min(test_mass_vector), max(test_mass_vector), length = 100)
  test_mean1 <- 70.00938
  test_mean2 <- 70.00962
  test_scale_factor1 <- 2
  test_scale_factor2 <- 4
  test_sigma1 <- 0.0001
  test_sigma2 <- 0.0002
  test_resol <- 140000
  
  expect_type(sum_curves(test_mean1, test_mean2, test_scale_factor1, test_scale_factor2, test_sigma1, test_sigma2, test_mass_vector_eq, test_mass_vector, test_resol), "list")
  expect_equal(sum_curves(test_mean1, test_mean2, test_scale_factor1, test_scale_factor2, test_sigma1, test_sigma2, test_mass_vector_eq, test_mass_vector, test_resol)$mean, 70.0095, tolerance = 0.0001, TRUE)
})

# test fit_optim
test_that("optimal peak fit can be found", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_int_vector <- c(9000, 16000, 19000, 15000, 6000)
  test_resol <- 140000
  # create output to test on:
  test_output <- list(mean = 70.0095, area = 5009.596, min = 70.00906, max = 70.00994)
  
  expect_type(fit_optim(test_mass_vector, test_int_vector, test_resol), "list")
  for (key in names(test_output)) {
    expect_equal(fit_optim(test_mass_vector, test_int_vector, test_resol)[[key]], test_output[[key]], tolerance = 0.0001, TRUE)
  }
})

# test optimize_1gaussian:
test_that("optimal value for m/z and scale_factor can be found", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_int_vector <- c(9000, 16000, 19000, 15000, 6000)
  test_sigma <- 0.0001
  test_query_mass <- median(test_mass_vector)
  test_scale_factor <- 2
  test_use_bounds <- FALSE
  
  expect_type(optimize_1gaussian(test_mass_vector, test_int_vector, test_sigma, test_query_mass, test_scale_factor, test_use_bounds), "list")
  expect_equal(optimize_1gaussian(test_mass_vector, test_int_vector, test_sigma, test_query_mass, test_scale_factor, test_use_bounds)$par[1], 70.00949, tolerance = 0.0001, TRUE)
  expect_equal(optimize_1gaussian(test_mass_vector, test_int_vector, test_sigma, test_query_mass, test_scale_factor, test_use_bounds)$par[2], 4.570024, tolerance = 0.0001, TRUE)
})

# test fit_1peak:
test_that("correct mean m/z and scale_factor can be found for 1 peak", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_int_vector <- c(9000, 16000, 19000, 15000, 6000)
  test_mass_vector_eq <- seq(min(test_mass_vector), max(test_mass_vector), length = 100)
  test_max_index <- which(test_int_vector == max(test_int_vector))
  test_resol <- 140000
  test_fit_quality <- 0.5
  test_use_bounds <- FALSE
  # create output to test on:
  test_output <- list(mean = 70.0095, scale_factor = 5.23334, sigma = 0.000126, qual = 0.24300)
  
  expect_type(fit_1peak(test_mass_vector_eq, test_mass_vector, test_int_vector, test_max_index, test_resol, test_fit_quality, test_use_bounds), "list")
  for (key in names(test_output)) {
    expect_equal(fit_1peak(test_mass_vector_eq, test_mass_vector, test_int_vector, test_max_index, test_resol, test_fit_quality, test_use_bounds)[[key]], test_output[[key]], tolerance = 0.0001, TRUE)
  }
})

# test optimize_2gaussians: 
test_that("optimal value for m/z and scale_factor can be found", {
  # create info for two peaks in a small region
  test_mass_vector_2peaks <- rep(70.00938, 14) + 1:14 * 0.00006
  test_int_vector_2peaks <- rep(c(2000, 9000, 16000, 19000, 15000, 6000, 3000), 2)
  test_mass_vector_eq <- seq(min(test_mass_vector_2peaks), max(test_mass_vector_2peaks), length = 100)
  test_max_index <- which(test_int_vector_2peaks == max(test_int_vector_2peaks))
  test_sigma1 <- 0.0001
  test_sigma2 <- 0.0002
  test_query_mass1 <- test_max_index[1]
  test_query_mass2 <- test_max_index[2]
  test_scale_factor1 <- 2
  test_scale_factor2 <- 3
  test_use_bounds <- FALSE
  
  expect_type(optimize_2gaussians(test_mass_vector_2peaks, test_int_vector_2peaks, test_sigma1, test_sigma2, test_query_mass1, test_scale_factor1, test_query_mass2, test_scale_factor2, test_use_bounds), "list")
  expect_equal(optimize_2gaussians(test_mass_vector_2peaks, test_int_vector_2peaks, test_sigma1, test_sigma2, test_query_mass1, test_scale_factor1, test_query_mass2, test_scale_factor2, test_use_bounds)$par[1], 4, tolerance = 0.1, TRUE)
})

# test fit_2peaks:
test_that("correct mean m/z and scale_factor can be found for 2 peaks", {
  # create info for two peaks in a small region
  test_mass_vector_2peaks <- rep(70.00938, 14) + 1:14 * 0.00006
  test_int_vector_2peaks <- rep(c(2000, 9000, 16000, 19000, 15000, 6000, 3000), 2)
  test_mass_vector_eq <- seq(min(test_mass_vector_2peaks), max(test_mass_vector_2peaks), length = 100)
  test_max_index <- which(test_int_vector_2peaks == max(test_int_vector_2peaks))
  test_resol <- 140000
  test_fit_quality <- 0.5
  test_use_bounds <- FALSE
  # create output to test on:
  test_output <- list(mean = c(70.00954, 70.00998), scale_factor = c(4.823598, 4.828072), sigma = c(0.0001257425, 0.0001257436), qual = 0.38341)
  
  expect_type(fit_2peaks(test_mass_vector_eq, test_mass_vector_2peaks, test_int_vector_2peaks, test_max_index, 
                         test_resol, test_use_bounds, test_fit_quality), "list")
  for (key in names(test_output)) {
    expect_equal(fit_2peaks(test_mass_vector_eq, test_mass_vector_2peaks, test_int_vector_2peaks, test_max_index, 
                            test_resol, test_use_bounds, test_fit_quality)[[key]], test_output[[key]], tolerance = 0.0001, TRUE)
  }
})

# test fit_gaussian(mass_vector2, mass_vector, int_vector, resol, force, use_bounds)
test_that("initial gaussian fit is done correctly", {
  # create peak info to test on:
  test_mass_vector <- rep(70.00938, 5) + 1:5 * 0.00006
  test_int_vector <- c(9000, 16000, 19000, 15000, 6000)
  test_mass_vector_eq <- seq(min(test_mass_vector), max(test_mass_vector), length = 100)
  test_resol <- 140000
  test_use_bounds <- FALSE
  test_force_nr <- 1
  
  expect_type(fit_gaussian(test_mass_vector_eq, test_mass_vector, test_int_vector, test_resol, test_force_nr, test_use_bounds), "list")
  expect_equal(fit_gaussian(test_mass_vector_eq, test_mass_vector, test_int_vector, test_resol, test_force_nr, test_use_bounds)$mean, 70.00956, tolerance = 0.00001, TRUE)
})

# test search_mzrange(ints_fullrange, resol, sample_name, peak_thresh)
test_that("all peak finding functions work together", {
  # enable snapshot
  local_edition(3)
  # create info for ten peaks separated by zero
  test_large_mass_vector <- rep(70.00938, 80) + 1:80 * 0.00006
  test_large_int_vector <- rep(c(2000, 9000, 16000, 19000, 15000, 6000, 3000, 0), 10)
  names(test_large_int_vector) <- test_large_mass_vector
  test_resol <- 140000
  test_sample_name <- "C1.1"
  test_peak_thresh <- 100
  
  expect_type(search_mzrange(test_large_int_vector, test_resol, test_sample_name, test_peak_thresh), "list")
  expect_snapshot(search_mzrange(test_large_int_vector, test_resol, test_sample_name, test_peak_thresh), error = FALSE)
})

# test fit_gaussian(mass_vector2, mass_vector, int_vector, resol, force, use_bounds)
test_that("fit gaussian")
