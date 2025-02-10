# functions for peak finding
search_mzrange <- function(ints_fullrange, resol, sample_name, peak_thresh) {
  #' Divide the full m/z range into regions of interest with min, max and mean m/z
  #'
  #' @param ints_fullrange: Named list of intensities (float)
  #' @param resol: Value for resolution (integer)
  #' @param sample_name: Sample name (string)
  #' @param peak_thresh: Value for noise level threshold (integer)
  #'
  #' @return allpeaks_values: list of m/z regions of interest

  # initialize list to store results for all peaks
  allpeaks_values <- list("mean" = NULL, "area" = NULL, "nr" = NULL,
                          "min" = NULL, "max" = NULL, "qual" = NULL, "spikes" = 0)
  
  # find indices where intensity is not equal to zero
  nonzero_positions <- as.vector(which(ints_fullrange != 0))

  # initialize
  # start position of the first peak
  start_index <- nonzero_positions[1]
  # maximum length of region of interest
  max_roi_length <- 15

  # find regions of interest
  for (running_index in 1:length(nonzero_positions)) {
    # find position of the end of a peak.
    if (running_index < length(nonzero_positions) && (nonzero_positions[running_index + 1] - nonzero_positions[running_index]) > 1) {
      end_index <- nonzero_positions[running_index]
      # get m/z values and intensities for this region of interest
      mass_vector <- as.numeric(names(ints_fullrange)[c(start_index:end_index)])
      int_vector <- as.vector(ints_fullrange[c(start_index:end_index)])
      # check if intensity is above threshold or the maximum intensity is NaN
      if (max(int_vector) < peak_thresh || is.nan(max(int_vector))) {
        # go to next region of interest
        start_index <- nonzero_positions[running_index + 1]
        next
      }
      # check if there are more intensities than maximum for region of interest
      if (length(int_vector) > max_roi_length) {
        print(length(int_vector))
        print(running_index)
        # trim lowest intensities to zero
        #int_vector[which(int_vector < min(int_vector) * 1.1)] <- 0
        # split the range into multiple sub ranges
        #sub_range <- int_vector
        #names(sub_range) <- mass_vector
        #allpeaks_values <- search_mzrange(sub_range, allpeaks_values, resol, 
        #                                  sample_name, peak_thresh)
      # A proper peak needs to have at least 3 intensities above threshold  
      } else if (length(int_vector) > 3) {
        # check if the sum of intensities is above zero. Why is this necessary?
        #if (sum(int_vector) == 0) next
       # define mass_diff as difference between last and first value of mass_vector
        # mass_diff <- mass_vector[length(mass_vector)] - mass_vector[1]
        # generate a second mass_vector with equally spaced m/z values
        mass_vector_eq <- seq(mass_vector[1], mass_vector[length(mass_vector)],
                            length = 10 * length(mass_vector))
        
        # Find the index in int_vector with the highest intensity
        # max_index <- which(int_vector == max(int_vector))
        # get initial fit values
        roi_values <- fit_gaussian(mass_vector_eq, mass_vector, int_vector, 
                                   resol, force_nr = length(max_index),
                                   use_bounds = FALSE)
        
        if (roi_values$qual[1] == 1) {
          # get optimized fit values
          roi_values <- fit_optim(mass_vector, int_vector, resol)
          # add region of interest to list of all peaks
          allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean)
          allpeaks_values$area <- c(allpeaks_values$area, roi_values$area)
          allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
          allpeaks_values$min <- c(allpeaks_values$min, roi_values$min)
          allpeaks_values$max <- c(allpeaks_values$max, roi_values$max)
          allpeaks_values$qual <- c(allpeaks_values$qual, 0)
          allpeaks_values$spikes <- allpeaks_values$spikes + 1
          
        } else {
          for (j in 1:length(roi_values$mean)){
            allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean[j])
            allpeaks_values$area <- c(allpeaks_values$area, roi_values$area[j])
            allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
            allpeaks_values$min <- c(allpeaks_values$min, roi_values$min[1])
            allpeaks_values$max <- c(allpeaks_values$max, roi_values$max[1])
            allpeaks_values$qual <- c(allpeaks_values$qual, roi_values$qual[1])
          }
        }
        
      } else {
        
        roi_values <- fit_optim(mass_vector, int_vector, resol)
        allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean)
        allpeaks_values$area <- c(allpeaks_values$area, roi_values$area)
        allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
        allpeaks_values$min <- c(allpeaks_values$min, roi_values$min)
        allpeaks_values$max <- c(allpeaks_values$max, roi_values$max)
        allpeaks_values$qual <- c(allpeaks_values$qual, 0)
        allpeaks_values$spikes <- allpeaks_values$spikes + 1
      }
    }
    start_index <- nonzero_positions[running_index + 1]
  }

  return(allpeaks_values)
}

fit_gaussian <- function(mass_vector_eq, mass_vector, int_vector,  
                         resol, force_nr, use_bounds) {
  #' Fit 1, 2, 3 or 4 Gaussian peaks in small region of m/z
  #'
  #' @param mass_vector_eq: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param resol: Value for resolution (integer)
  #' @param force_nr: Number of local maxima in int_vector (integer)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  # Find the index in int_vector with the highest intensity
  max_index <- which(int_vector == max(int_vector))
  
  # Initialise
  peak_mean <- NULL
  peak_area <- NULL
  peak_qual <- NULL
  peak_min <- NULL
  peak_max <- NULL
  fit_quality1 <- 0.15
  fit_quality <- 0.2
  # set initial value for scale factor
  scale_factor <- 2
  
  # One local maximum:
  if (force_nr == 1) {
    # determine fit values for 1 Gaussian peak 
    fit_values <- fit_1peak(mass_vector_eq, mass_vector, int_vector, max_index, resol,
                            fit_quality1, use_bounds)
    
    # test if the mean is outside the m/z range
    if (fit_values$mean[1] < mass_vector[1] || fit_values$mean[1] > mass_vector[length(mass_vector)]) {
      # run this function again with fixed boundaries
      return(fit_gaussian(mass_vector_eq, mass_vector, int_vector, resol,
                          force_nr = 1, use_bounds = TRUE))
    } else {
      # test if the fit is bad
      if (fit_values$qual > fit_quality1) {
        # Try to fit two curves; find two local maxima. NB: max_index (now new_index) removed from fit_gaussian
        new_index <- which(diff(sign(diff(int_vector))) == -2) + 1
        # test if there are two indices in new_index
        if (length(new_index) != 2) {
          new_index <- round(length(mass_vector) / 3)
          new_index <- c(new_index, 2 * new_index)
        }
        # run this function again with two local maxima
        return(fit_gaussian(mass_vector_eq, mass_vector, int_vector, 
                            resol, force_nr = 2, use_bounds = FALSE))
        # good fit
      } else {
        peak_mean <- c(peak_mean, fit_values$mean)
        peak_area <- c(peak_area, estimate_area(fit_values$mean, resol, fit_values$scale_factor,
                                                fit_values$sigma))
        peak_qual <- fit_values$qual
        peak_min <- mass_vector[1]
        peak_max <- mass_vector[length(mass_vector)]
      }
    }
    
    #### Two local maxima; need at least 6 data points for this ####
  } else if (force_nr == 2 && (length(mass_vector) > 6)) {
    # determine fit values for 2 Gaussian peaks
    fit_values <- fit_2peaks(mass_vector_eq, mass_vector, int_vector, max_index, resol,
                             use_bounds, fit_quality)
    # test if one of the means is outside the m/z range
    if (fit_values$mean[1] < mass_vector[1] || fit_values$mean[1] > mass_vector[length(mass_vector)] ||
        fit_values$mean[2] < mass_vector[1] || fit_values$mean[2] > mass_vector[length(mass_vector)]) {
      # check if fit quality is bad
      if (fit_values$qual > fit_quality) {
        # run this function again with fixed boundaries
        return(fit_gaussian(mass_vector_eq, mass_vector, int_vector, resol,
                            force_nr = 2, use_bounds = TRUE))
      } else {
        # check which mean is outside range and remove it from the list of means
        # NB: peak_mean and other variables have not been given values from 2-peak fit yet!
        for (i in 1:length(fit_values$mean)){
          if (fit_values$mean[i] < mass_vector[1] || fit_values$mean[i] > mass_vector[length(mass_vector)]) {
            peak_mean <- c(peak_mean, -i)
            peak_area <- c(peak_area, -i)
          } else {
            peak_mean <- c(peak_mean, fit_values$mean[i])
            peak_area <- c(peak_area, fit_values$area[i])
          }
        }
        peak_qual <- fit_values$qual
        peak_min <- mass_vector[1]
        peak_max <- mass_vector[length(mass_vector)]
      }
      # if all means are within range
    } else {
      # check for bad fit
      if (fit_values$qual > fit_quality) {
        # Try to fit three curves; find three local maxima
        new_index <- which(diff(sign(diff(int_vector))) == -2) + 1
        # test if there are three indices in new_index
        if (length(new_index) != 3) {
          new_index <- round(length(mass_vector) / 4)
          new_index <- c(new_index, 2 * new_index, 3 * new_index)
        }
        # run this function again with three local maxima
        return(fit_gaussian(mass_vector_eq, mass_vector, int_vector, 
                            resol, force_nr = 3, use_bounds = FALSE))
      # good fit, all means are within m/z range
      } else {
        # check if means are within 3 ppm and sum if so
        tmp <- fit_values$qual
        nr_means_new <- -1
        nr_means <- length(fit_values$mean)
        while (nr_means != nr_means_new) {
          nr_means <- length(fit_values$mean)
          fit_values <- within_ppm(fit_values$mean, fit_values$scale_factor, fit_values$sigma, fit_values$area, 
                                   mass_vector_eq, mass_vector, ppm = 4, resol)
          nr_means_new <- length(fit_values$mean)
        }
        # restore original quality score
        fit_values$qual <- tmp
        
        for (i in 1:length(fit_values$mean)){
          peak_mean <- c(peak_mean, fit_values$mean[i])
          peak_area <- c(peak_area, fit_values$area[i])
        }
        peak_qual <- fit_values$qual
        peak_min <- mass_vector[1]
        peak_max <- mass_vector[length(mass_vector)]
      }
    }
    
  } else { # More than two local maxima; fit 1 peak.
    fit_quality1 <- 0.40
    use_bounds <- TRUE
    max_index <- which(int_vector == max(int_vector))
    fit_values <- fit_1peak(mass_vector_eq, mass_vector, int_vector, max_index, resol,
                            fit_quality1, use_bounds)
    # check for bad fit
    if (fit_values$qual > fit_quality1) {
      # get fit values from fit_optim
      fit_values <- fit_optim(mass_vector, int_vector, resol)
      peak_mean <- c(peak_mean, fit_values$mean)
      peak_area <- c(peak_area, fit_values$area)
      peak_min <- fit_values$min
      peak_max <- fit_values$max
      peak_qual <- 0
    } else {
      peak_mean <- c(peak_mean, fit_values$mean)
      peak_area <- c(peak_area, estimate_area(fit_values$mean, resol, fit_values$scale_factor, fit_values$sigma))
      peak_qual <- fit_values$qual
      peak_min <- mass_vector[1]
      peak_max <- mass_vector[length(mass_vector)]
    }
  }
  
  # put all values for this region of interest into a list
  roi_value_list <- list("mean" = peak_mean,
                         "area" = peak_area,
                         "qual" = peak_qual,
                         "min" = peak_min,
                         "max" = peak_max)
  return(roi_value_list)
}


fit_1peak <- function(mass_vector_eq, mass_vector, int_vector, max_index, 
                      resol, fit_quality, use_bounds) {
  #' Fit 1 Gaussian peak in small region of m/z
  #'
  #' @param mass_vector_eq: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param max_index: Index in int_vector with the highest intensity (integer)
  #' @param resol: Value for resolution (integer)
  #' @param fit_quality: Value indicating quality of fit of Gaussian curve (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  # set initial value for scale_factor
  scale_factor <- 2
  
  if (length(int_vector) < 3) {
    message("Range too small, no fit possible!")
  } else {
    if (length(int_vector) == 4) {
      # fit 1 peak
      weighted_mu <- weighted.mean(mass_vector, int_vector)
      sigma <- get_stdev(mass_vector, int_vector)
      fitted_peak <- optimize_1gaussian(mass_vector, int_vector, sigma, weighted_mu, scale_factor, use_bounds)
    } else {
      # set range vector
      if ((length(mass_vector) - length(max_index)) < 2) {
        range1 <- c((length(mass_vector) - 4) : length(mass_vector))
      } else if (length(max_index) < 2) {
        range1 <- c(1:5)
      } else {
        range1 <- seq(from = (max_index[1] - 2), to = (max_index[1] + 2))
      }
      # remove zero at the beginning of range1
      if (range1[1] == 0) {
        range1 <- range1[-1]
      }
      # remove NA
      if (length(which(is.na(int_vector[range1]))) != 0) {
        range1 <- range1[-which(is.na(int_vector[range1]))]
      }
      # fit 1 peak
      weighted_mu <- weighted.mean(mass_vector[range1], int_vector[range1])
      sigma <- get_stdev(mass_vector[range1], int_vector[range1])
      fitted_peak <- optimize_1gaussian(mass_vector, int_vector, sigma, weighted_mu, scale_factor, use_bounds)
    }
    
    fitted_mz <- fitted_peak$par[1]
    fitted_nr <- fitted_peak$par[2]
    
    # get new value for fit quality and scale_factor
    fq_new <- get_fit_quality(mass_vector, int_vector, fitted_mz, resol, fitted_nr, sigma)$fq_new
    scale_factor_new <- 1.2 * scale_factor
    
    # bad fit
    if (fq_new > fit_quality) {
      # optimize scaling factor
      fq <- 0
      scale_factor <- 0
      if (sum(int_vector) > sum(params1[2] * dnorm(mass_vector, fitted_mz, sigma))) {
        # increase scale_factor until convergence
        while ((round(fq, digits = 3) != round(fq_new, digits = 3)) && (scale_factor_new < 10000)) {
          fq <- fq_new
          scale_factor <- scale_factor_new
          # fit 1 peak
          fitted_peak <- optimize_1gaussian(mass_vector, int_vector, sigma, weighted_mu, scale_factor, use_bounds)
          params1 <- fitted_peak$par
          # get new value for fit quality and scale_factor
          fq_new <- get_fit_quality(mass_vector, int_vector, fitted_mz, resol, fitted_nr, sigma)$fq_new
          scale_factor_new <- 1.2 * scale_factor
        }
      } else {
        # decrease scale_factor until convergence
        while ((round(fq, digits = 3) != round(fq_new, digits = 3)) && (scale_factor_new < 10000)) {
          fq <- fq_new
          scale_factor <- scale_factor_new
          # fit 1 peak
          fitted_peak <- optimize_1gaussian(mass_vector, int_vector, sigma, weighted_mu, scale_factor, use_bounds)
          params1 <- fitted_peak$par
          # get new value for fit quality and scale_factor
          fq_new <- get_fit_quality(mass_vector, int_vector, fitted_mz, resol, fitted_nr, sigma)$fq_new
          scale_factor_new <- 0.8 * scale_factor
          print(scale_factor_new)
        }
      }
      # use optimized scale_factor factor to fit 1 peak
      if (fq < fq_new) {
        fitted_peak <- optimize_1gaussian(mass_vector, int_vector, sigma, weighted_mu, scale_factor, use_bounds)
        fitted_mz <- fitted_peak$par[1]
        fitted_nr <- fitted_peak$par[2]
        fq_new <- fq
      }
    }
  }
  
  roi_value_list <- list("mean" = fitted_mz, "scale_factor" = fitted_nr, "sigma" = sigma, "qual" = fq_new)
  return(roi_value_list)
}

fit_2peaks <- function(mass_vector_eq, mass_vector, int_vector, max_index, resol, use_bounds = FALSE,
                       fit_quality) {
  #' Fit 2 Gaussian peaks in a small region of m/z
  #'
  #' @param mass_vector_eq: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param max_index: Index in int_vector with the highest intensity (integer)
  #' @param resol: Value for resolution (integer)
  #' @param fit_quality: Value indicating quality of fit of Gaussian curve (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  peak_mean <- NULL
  peak_area <- NULL
  peak_scale <- NULL
  peak_sigma <- NULL
  
  # set range vectors for first peak
  range1 <- seq(from = (max_index[1] - 2), to = (max_index[1] + 2))
  # remove zero at the beginning of range1
  if (range1[1] == 0) {
    range1 <- range1[-1]
  }
  # set range vectors for second peak
  range2 <- seq(from = (max_index[2] - 2), to = (max_index[2] + 2))
  # if range2 ends outside mass_vector, shorten it
  if (length(mass_vector) < range2[length(range2)]) {
    range2 <- range2[-length(range2)]
  }
  # check whether the two ranges overlap
  range1 <- check_overlap(range1, range2)[[1]]
  range2 <- check_overlap(range1, range2)[[2]]
  # check for negative or 0 values in range1 or range2
  remove <- which(range1 < 1)
  if (length(remove) > 0) {
    range1 <- range1[-remove]
  }
  remove <- which(range2 < 1)
  if (length(remove) > 0) {
    range2 <- range2[-remove]
  }
  # remove NA
  if (length(which(is.na(int_vector[range1]))) != 0) {
    range1 <- range1[-which(is.na(int_vector[range1]))]
  }
  if (length(which(is.na(int_vector[range2]))) != 0) {
    range2 <- range2[-which(is.na(int_vector[range2]))]
  }
  
  # fit 2 peaks, first separately, then together
  weighted_mu1 <- weighted.mean(mass_vector[range1], int_vector[range1])
  sigma1 <- get_stdev(mass_vector[range1], int_vector[range1])
  fitted_peak1 <- optimize_1gaussian(mass_vector[range1], int_vector[range1], sigma1, weighted_mu1, scale_factor, use_bounds)
  fitted_mz1 <- fitted_peak1$par[1]
  fitted_nr1 <- fitted_peak1$par[2]
  # second peak
  weighted_mu2 <- weighted.mean(mass_vector[range2], int_vector[range2])
  sigma2 <- get_stdev(mass_vector[range2], int_vector[range2])
  fitted_peak2 <- optimize_1gaussian(mass_vector[range2], int_vector[range2], sigma2, weighted_mu2, scale_factor, use_bounds)
  fitted_mz2 <- fitted_peak2$par[1]
  fitted_nr2 <- fitted_peak2$par[2]
  # combined
  fitted_2peaks <- optimize_2gaussians(mass_vector, int_vector, sigma1, sigma2, 
                                       fitted_mz1, fitted_nr1, 
                                       fitted_mz2, fitted_nr2, use_bounds)
  fitted_2peaks_mz1 <- fitted_2peaks$par[1]
  fitted_2peaks_nr1 <- fitted_2peaks$par[2]
  fitted_2peaks_mz2 <- fitted_2peaks$par[3]
  fitted_2peaks_nr2 <- fitted_2peaks$par[4]
  
  # get fit quality
  if (is.null(sigma2)) {
    sigma2 <- sigma1
  }
  sum_fit <- (fitted_2peaks_nr1 * dnorm(mass_vector, fitted_2peaks_mz1, sigma1)) +
    (fitted_2peaks_nr2 * dnorm(mass_vector, fitted_2peaks_mz2, sigma2))
  lowest_mz <- sort(c(fitted_2peaks_mz1, fitted_2peaks_mz2))[1]
  fq_new <- get_fit_quality(mass_vector, int_vector, lowest_mz,
                        resol, sum_fit = sum_fit)$fq_new
  
  # get parameter values
  area1 <- estimate_area(fitted_2peaks_mz1, resol, fitted_2peaks_nr1, sigma1)
  area2 <- estimate_area(fitted_2peaks_mz2, resol, fitted_2peaks_nr2, sigma2)
  peak_area <- c(peak_area, area1, area2)
  peak_mean <- c(peak_mean, fitted_2peaks_mz1, fitted_2peaks_mz2)
  peak_scale <- c(peak_scale, fitted_2peaks_nr1, fitted_2peaks_nr2)
  peak_sigma <- c(peak_sigma, sigma1, sigma2)
  
  roi_value_list <- list("mean" = peak_mean, "scale_factor" = peak_scale, "sigma" = peak_sigma, "area" = peak_area, "qual" = fq_new)
  return(roi_value_list)
}

optimize_1gaussian <- function(mass_vector, int_vector, sigma, query_mass, scale_factor, use_bounds) {
  #' Fit a Gaussian curve for a peak with given parameters
  #'
  #' @param mass_vector: Vector of masses (float)
  #' @param int_vector: Vector of intensities (float)
  #' @param sigma: Value for width of the peak (float)
  #' @param query_mass: Value for mass at center of peak (float)
  #' @param scale_factor: Value for scaling intensities (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return opt_fit: list of parameters and values describing the optimal fit
  
  # define optimization function for optim based on normal distribution
  opt_f <- function(params) {
    d <- params[2] * dnorm(mass_vector, mean = params[1], sd = sigma)
    sum((d - int_vector) ^ 2)
  }
  if (use_bounds) {
    # determine lower and upper boundaries
    lower <- c(mass_vector[1], 0, mass_vector[1], 0)
    upper <- c(mass_vector[length(mass_vector)], Inf, mass_vector[length(mass_vector)], Inf)
    # get optimal value for fitted Gaussian curve
    opt_fit <- optim(c(as.numeric(query_mass), as.numeric(scale_factor)),
                     opt_f, control = list(maxit = 10000), method = "L-BFGS-B",
                     lower = lower, upper = upper)
  } else {
    opt_fit <- optim(c(as.numeric(query_mass), as.numeric(scale_factor)),
                     opt_f, control = list(maxit = 10000))
  }
  return(opt_fit)
}

optimize_2gaussians <- function(mass_vector, int_vector, sigma1, sigma2,
                           query_mass1, scale_factor1,
                           query_mass2, scale_factor2, use_bounds) {
  #' Fit two Gaussian curves for a peak with given parameters
  #'
  #' @param mass_vector: Vector of masses (float)
  #' @param int_vector: Vector of intensities (float)
  #' @param sigma1: Value for width of the first peak (float)
  #' @param sigma2: Value for width of the second peak (float)
  #' @param query_mass1: Value for mass at center of first peak (float)
  #' @param scale_factor1: Value for scaling intensities for first peak (float)
  #' @param query_mass2: Value for mass at center of second peak (float)
  #' @param scale_factor2: Value for scaling intensities for second peak (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return opt_fit: list of parameters and values describing the optimal fit
  
  # define optimization function for optim based on normal distribution
  opt_f <- function(params) {
    d <- params[2] * dnorm(mass_vector, mean = params[1], sd = sigma1) +
      params[4] * dnorm(mass_vector, mean = params[3], sd = sigma2)
    sum((d - int_vector) ^ 2)
  }
  
  if (use_bounds) {
    # determine lower and upper boundaries
    lower <- c(mass_vector[1], 0, mass_vector[1], 0)
    upper <- c(mass_vector[length(mass_vector)], Inf, mass_vector[length(mass_vector)], Inf)
    # get optimal value for 2 fitted Gaussian curves
    if (is.null(query_mass2) && is.null(scale_factor2) && is.null(sigma2)) {
      sigma2 <- sigma1
      opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale_factor1),
                         as.numeric(query_mass1), as.numeric(scale_factor1)),
                       opt_f, control = list(maxit = 10000),
                       method = "L-BFGS-B", lower = lower, upper = upper)
    } else {
      opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale_factor1),
                         as.numeric(query_mass2), as.numeric(scale_factor2)),
                       opt_f, control = list(maxit = 10000),
                       method = "L-BFGS-B", lower = lower, upper = upper)
    }
  } else {
    if (is.null(query_mass2) && is.null(scale_factor2) && is.null(sigma2)) {
      sigma2 <- sigma1
      opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale_factor1),
                         as.numeric(query_mass1), as.numeric(scale_factor1)),
                       opt_f, control = list(maxit = 10000))
    } else {
      opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale_factor1),
                         as.numeric(query_mass2), as.numeric(scale_factor2)),
                       opt_f, control = list(maxit = 10000))
    }
  }
  return(opt_fit)
}

get_stdev <- function(mass_vector, int_vector, resol = 140000) {
  #' Calculate standard deviation to determine width of a peak
  #'
  #' @param mass_vector: Vector of 3 mass values (float)
  #' @param int_vector: Vector of 3 intensities (float)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return stdev: Value for standard deviation
  # find maximum intensity in vector
  max_index <- which(int_vector == max(int_vector))
  # find corresponding mass at maximum intensity
  max_mass <- mass_vector[max_index]
  # calculate resolution at given m/z value
  resol_mz <- resol * (1 / sqrt(2) ^ (log2(max_mass / 200)))
  # calculate full width at half maximum
  fwhm <- max_mass / resol_mz
  # calculate standard deviation
  stdev <- (fwhm / 2) * 0.85
  return(stdev)
}

get_fit_quality <- function(mass_vector, int_vector, mu_first, resol, scale_factor = NULL, sigma = NULL, sum_fit = NULL) {
  #' Get fit quality for 1 Gaussian peak in small region of m/z
  #'
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param mu_first: Value for first peak (float)
  #' @param scale_factor: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param sum_fit: Value indicating quality of fit of Gaussian curve (float)
  #'
  #' @return list_params: list of parameters indicating quality of fit (list)
  if (is.null(sum_fit)) {
    mass_vector_int <- mass_vector
    int_vector_int <- int_vector
    # get new fit quality
    fq_new <- mean(abs((scale_factor * dnorm(mass_vector_int, mu_first, sigma)) - int_vector_int) / 
                     rep((max(scale_factor * dnorm(mass_vector_int, mu_first, sigma)) / 2), length(mass_vector_int)))
  } else {
    sum_fit_int <- sum_fit
    int_vector_int <- int_vector
    mass_vector_int <- mass_vector
    # get new fit quality
    fq_new <- mean(abs(sum_fit_int - int_vector_int) / rep(max(sum_fit_int) /2, length(sum_fit_int)))
  }
  
  # Prevent division by 0
  if (is.nan(fq_new)) fq_new <- 1 
  
  list_params <- list("fq_new" = fq_new, "x_int" = mass_vector_int, "y_int" = int_vector_int)
  return(list_params)
}

estimate_area <- function(mass_max, resol, scale_factor, sigma) {
  #' Estimate area of Gaussian curve
  #'
  #' @param mass_max: Value for m/z at maximum intensity of a peak (float)
  #' @param resol: Value for resolution (integer)
  #' @param scale_factor: Value for peak width (float)
  #' @param sigma: Value for standard deviation (float)
  #'
  #' @return area_curve: Value for area under the Gaussian curve (float)
  
  # calculate width of peak at half maximum
  fwhm <- get_fwhm(mass_max, resol)
  
  # generate a mass_vector with equally spaced m/z values
  mz_min <- mass_max - 2 * fwhm
  mz_max <- mass_max + 2 * fwhm
  mz_range <- mz_max - mz_min 
  mass_vector_eq <- seq(mz_min, mz_max, length = 100)
  
  # estimate area under the curve
  area_curve <- sum(scale_factor * dnorm(mass_vector_eq, mass_max, sigma)) / 100
  
  return(area_curve)  
}

get_fwhm <- function(query_mass, resol) {
  #' Calculate fwhm (full width at half maximum intensity) for a peak
  #'
  #' @param query_mass: Value for mass (float)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return fwhm: Value for full width at half maximum (float)

  # set aberrant values of query_mass to default value of 200
  if (is.na(query_mass)) {
    query_mass <- 200
  } 
  if (query_mass < 0) {
    query_mass <- 200
  }
  # calculate resolution at given m/z value
  resol_mz <- resol * (1 / sqrt(2) ^ (log2(query_mass / 200)))
  # calculate full width at half maximum
  fwhm <- query_mass / resol_mz
  
  return(fwhm)
}

fit_optim <- function(mass_vector, int_vector, resol) {
  #' Determine optimized fit of Gaussian curve to small region of m/z
  #'
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Vector of intensities for a region of interest (float)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  # initial value for scale_factor
  scale_factor <- 1.5

  # Find the index in int_vector with the highest intensity
  max_index <- which(int_vector == max(int_vector))[1]
  mass_max <- mass_vector[max_index]
  int_max <- int_vector[max_index]
  # get peak width
  fwhm <- get_fwhm(mass_max, resol)
  # simplify the peak shape: represent it by a triangle
  mass_max_simple <- c(mass_max - scale_factor * fwhm, mass_max, mass_max + scale_factor * fwhm)
  int_max_simple <- c(0, int_max, 0)
  
  # define mass_diff as difference between last and first value of mass_max_simple
  # mass_diff <- mass_max_simple[length(mass_max_simple)] - mass_max_simple[1]
  # generate a second mass_vector with equally spaced m/z values
  mass_vector_eq <- seq(mass_max_simple[1], mass_max_simple[length(mass_max_simple)], 
                      length = 100 * length(mass_max_simple))
  sigma <- get_stdev(mass_vector_eq, int_max_simple)
  # define optimization function for optim based on normal distribution
  opt_f <- function(p, mass_vector, int_vector, sigma, mass_max) {
    curve <- p * dnorm(mass_vector, mass_max, sigma)
    return((max(curve) - max(int_vector))^2)
  }
  
  # get optimal value for fitted Gaussian curve
  opt_fit <- optimize(opt_f, c(0, 100000), tol = 0.0001, mass_vector, int_vector, sigma, mass_max)
  scale_factor <- opt_fit$minimum
  
  # get an estimate of the area under the peak
  area <- estimate_area(mass_max, resol, scale_factor, sigma)
  # put all values for this region of interest into a list
  roi_value_list <- list("mean" = mass_max,
                         "area" = area,
                         "min" = mass_vector_eq[1],
                         "max" = mass_vector_eq[length(mass_vector_eq)])
  return(roi_value_list)
}

within_ppm <- function(mean, scale_factor, sigma, area, mass_vector_eq, mass_vector, ppm = 4, resol) {
  #' Test whether two mass ranges are within ppm distance of each other
  #'
  #' @param mean: Value for mean m/z (float)
  #' @param scale_factor: Initial value used to estimate scaling parameter (integer)
  #' @param sigma: Value for standard deviation (float)
  #' @param area: Value for area under the curve (float)
  #' @param mass_vector_eq: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param ppm: Value for distance between two values of mass (integer)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return list_params: list of parameters indicating quality of fit (list)
  
  # sort
  index <- order(mean)
  mean <- mean[index]
  scale_factor <- scale_factor[index]
  sigma <- sigma[index]
  area <- area[index]
  
  summed <- NULL
  remove <- NULL
  
  if (length(mean) > 1) {
    for (i in 2:length(mean)) {
      if ((abs(mean[i - 1] - mean[i]) / mean[i - 1]) * 10^6 < ppm) {
        
        # avoid double occurrance in sum
        if ((i - 1) %in% summed) next
        
        result_values <- sum_curves(mean[i - 1], mean[i], scale_factor[i - 1], scale_factor[i], sigma[i - 1], sigma[i],
                                    mass_vector_eq, mass_vector, resol)
        summed <- c(summed, i - 1, i)
        if (is.nan(result_values$mean)) result_values$mean <- 0
        mean[i - 1] <- result_values$mean
        mean[i] <- result_values$mean
        area[i - 1] <- result_values$area
        area[i] <- result_values$area
        scale_factor[i - 1] <- result_values$scale_factor
        scale_factor[i] <- result_values$scale_factor
        sigma[i - 1] <- result_values$sigma
        sigma[i] <- result_values$sigma
        
        remove <- c(remove, i)
      }
    }
  }
  
  if (length(remove) != 0) {
    mean <- mean[-c(remove)]
    area <- area[-c(remove)]
    scale_factor <- scale_factor[-c(remove)]
    sigma <- sigma[-c(remove)]
  }
  
  list_params <- list("mean" = mean, "area" = area, "scale_factor" = scale_factor, "sigma" = sigma, "qual" = NULL)
  return(list_params)
}

sum_curves <- function(mean1, mean2, scale_factor1, scale_factor2, sigma1, sigma2, mass_vector_eq, mass_vector, resol) {
  #' Sum two curves
  #'
  #' @param mean1: Value for mean m/z of first peak (float)
  #' @param mean2: Value for mean m/z of second peak (float)
  #' @param scale_factor1: Initial value used to estimate scaling parameter for first peak (integer)
  #' @param scale_factor2: Initial value used to estimate scaling parameter for second peak (integer)
  #' @param sigma1: Value for standard deviation for first peak (float)
  #' @param sigma2: Value for standard deviation for second peak (float)
  #' @param mass_vector_eq: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return list_params: list of parameters indicating quality of fit (list)
  
  sum_fit <- (scale_factor1 * dnorm(mass_vector_eq, mean1, sigma1)) + (scale_factor2 * dnorm(mass_vector_eq, mean2, sigma2))
  
  mean1_plus2 <- weighted.mean(c(mean1, mean2), c(max(scale_factor1 * dnorm(mass_vector_eq, mean1, sigma1)),
                                                  max(scale_factor2 * dnorm(mass_vector_eq, mean2, sigma2))))
  
  # get new values for parameters
  fwhm <- get_fwhm(mean1_plus2, resol)
  area <- max(sum_fit)
  scale_factor <- scale_factor1 + scale_factor2
  sigma <- (fwhm / 2) * 0.85
  
  list_params <- list("mean" = mean1_plus2, "area" = area, "scale_factor" = scale_factor, "sigma" = sigma)
  return(list_params)
}

check_overlap <- function(range1, range2) {
  #' Modify range1 and range2 in case of overlap
  #'
  #' @param range1: Vector of m/z values for first peak (float)
  #' @param range2: Vector of m/z values for second peak (float)
  #'
  #' @return new_ranges: list of two ranges (list)
  
  # Check for overlap
  if (length(intersect(range1, range2)) == 2) {
    if (length(range1) >= length(range2)) {
      range1 <- range1[-length(range1)]  
    } else {
      range2 <- range2[-1]
    }
  } else if (length(intersect(range1, range2)) == 3) {
    range1 <- range1[-length(range1)]  
    range2 <- range2[-1]
  }
  new_ranges <- list("range1" = range1, "range2" = range2)
  return(new_ranges)
}

