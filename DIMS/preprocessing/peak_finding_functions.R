# functions for peak finding
search_mzrange <- function(ints_fullrange, resol, sample_name, scanmode, peak_thresh) {
  #' Divide the full m/z range into regions of interest with min, max and mean m/z
  #'
  #' @param ints_fullrange: Named list of intensities (float)
  #' @param resol: Value for resolution (integer)
  #' @param sample_name: Sample name (string)
  #' @param scanmode: Scan mode, positive or negative (string)
  #' @param peak_thresh: Value for noise level threshold (integer)
  #'
  #' @return allpeaks_values: list of m/z regions of interest

  # Number used to calculate area under Gaussian curve
  int_factor <- 1 * 10^5 
  
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
        #                                  sample_name, scanmode, peak_thresh)
      # A proper peak needs to have at least 3 intensities above threshold  
      } else if (length(int_vector) > 3) {
        # check if the sum of intensities is above zero. Why is this necessary?
        #if (sum(int_vector) == 0) next
       # define mass_diff as difference between last and first value of mass_vector
        mass_diff <- mass_vector[length(mass_vector)] - mass_vector[1]
        # generate a second mass_vector with equally spaced m/z values
        mass_vector2 <- seq(mass_vector[1], mass_vector[length(mass_vector)],
                            length = mass_diff * int_factor)
        
        # Find the index in int_vector with the highest intensity
        # max_index <- which(int_vector == max(int_vector))
        # get initial fit values
        roi_values <- fit_gaussian(mass_vector2, mass_vector, int_vector, # max_index,
                                   resol, force = length(max_index),
                                   use_bounds = FALSE, scanmode)
        
        if (roi_values$qual[1] == 1) {
          # get optimized fit values
          roi_values <- fit_optim(mass_vector, int_vector, resol, scanmode)
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
        
        roi_values <- fit_optim(mass_vector, int_vector, resol, scanmode)
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

# remove plot sections (commented out)
fit_gaussian <- function(mass_vector2, mass_vector, int_vector,  
                         resol, force, use_bounds, scanmode) {
  #' Fit 1, 2, 3 or 4 Gaussian peaks in small region of m/z
  #'
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param resol: Value for resolution (integer)
  #' @param force: Number of local maxima in int_vector (integer)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #' @param scanmode: Scan mode, positive or negative (string)
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  # Number used to calculate area under Gaussian curve
  int_factor <- 1 * 10^5 
  
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
  
  # One local maximum:
  if (force == 1) {
    # determine fit values for 1 Gaussian peak (mean, scale, sigma, qual)
    fit_values <- fit_1peak(mass_vector2, mass_vector, int_vector, max_index, resol,
                            fit_quality1, use_bounds)
    
    # set initial value for scale factor
    scale <- 2
    # test if the mean is outside the m/z range
    if (fit_values$mean[1] < mass_vector[1] || fit_values$mean[1] > mass_vector[length(mass_vector)]) {
      # run this function again with fixed boundaries
      return(fit_gaussian(mass_vector2, mass_vector, int_vector, resol,
                          force = 1, use_bounds = TRUE, scanmode))
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
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, 
                            resol, force = 2, use_bounds = FALSE, scanmode))
        # good fit
      } else {
        peak_mean <- c(peak_mean, fit_values$mean)
        peak_area <- c(peak_area, estimate_area(fit_values$mean, resol, fit_values$scale,
                                                fit_values$sigma))
        peak_qual <- fit_values$qual
        peak_min <- mass_vector[1]
        peak_max <- mass_vector[length(mass_vector)]
      }
    }
    
    #### Two local maxima; need at least 6 data points for this ####
  } else if (force == 2 && (length(mass_vector) > 6)) {
    # determine fit values for 2 Gaussian peaks (mean, scale, sigma, qual)
    fit_values <- fit_2peaks(mass_vector2, mass_vector, int_vector, max_index, resol,
                             use_bounds, fit_quality)
    # test if one of the means is outside the m/z range
    if (fit_values$mean[1] < mass_vector[1] || fit_values$mean[1] > mass_vector[length(mass_vector)] ||
        fit_values$mean[2] < mass_vector[1] || fit_values$mean[2] > mass_vector[length(mass_vector)]) {
      # check if fit quality is bad
      if (fit_values$qual > fit_quality) {
        # run this function again with fixed boundaries
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, resol,
                            force = 2, use_bounds = TRUE, scanmode))
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
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, 
                            resol, force = 3, use_bounds = FALSE, scanmode))
        # good fit, all means are within m/z range
      } else {
        # check if means are within 3 ppm and sum if so
        tmp <- fit_values$qual
        nr_means_new <- -1
        nr_means <- length(fit_values$mean)
        while (nr_means != nr_means_new) {
          nr_means <- length(fit_values$mean)
          fit_values <- within_ppm(fit_values$mean, fit_values$scale, fit_values$sigma, fit_values$area, 
                                   mass_vector2, mass_vector, ppm = 4, resol)
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
    scale <- 2
    fit_quality1 <- 0.40
    use_bounds <- TRUE
    max_index <- which(int_vector == max(int_vector))
    fit_values <- fit_1peak(mass_vector2, mass_vector, int_vector, max_index, resol,
                            fit_quality1, use_bounds)
    # check for bad fit
    if (fit_values$qual > fit_quality1) {
      # get fit values from fit_optim
      fit_values <- fit_optim(mass_vector, int_vector, resol, scanmode)
      peak_mean <- c(peak_mean, fit_values$mean)
      peak_area <- c(peak_area, fit_values$area)
      peak_min <- fit_values$min
      peak_max <- fit_values$max
      peak_qual <- 0
    } else {
      peak_mean <- c(peak_mean, fit_values$mean)
      peak_area <- c(peak_area, estimate_area(fit_values$mean, resol, fit_values$scale, fit_values$sigma))
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


fit_1peak <- function(mass_vector2, mass_vector, int_vector, max_index, 
                      resol, fit_quality, use_bounds) {
  #' Fit 1 Gaussian peak in small region of m/z
  #'
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param max_index: Index in int_vector with the highest intensity (integer)
  #' @param resol: Value for resolution (integer)
  #' @param fit_quality: Value indicating quality of fit of Gaussian curve (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  # set initial value for scale
  scale <- 2
  
  if (length(int_vector) < 3) {
    message("Range too small, no fit possible!")
  } else {
    if ((length(int_vector) == 4)) {
      # fit 1 peak
      mu <- weighted.mean(mass_vector, int_vector)
      sigma <- get_stdev(mass_vector, int_vector)
      fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, use_bounds)
    } else {
      # set range vector
      if ((length(mass_vector) - length(max_index)) < 2) {
        range1 <- c((length(mass_vector) - 4) : length(mass_vector))
      } else if (length(max_index) < 2) {
        range1 <- c(1:5)
      } else {
        range1 <- c(max_index[1] - 2, max_index[1] - 1, max_index[1], max_index[1] + 1, max_index[1] + 2)
      }
      if (range1[1] == 0) range1 <- range1[-1]
      # remove NA
      if (length(which(is.na(int_vector[range1]))) != 0) {
        range1 <- range1[-which(is.na(int_vector[range1]))]
      }
      # fit 1 peak
      mu <- weighted.mean(mass_vector[range1], int_vector[range1])
      sigma <- get_stdev(mass_vector[range1], int_vector[range1])
      fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, use_bounds)
    }
    
    p1 <- fitted_peak$par
    
    # get new value for fit quality and scale
    fq_new <- get_fit_quality(mass_vector, int_vector, p1[1], p1[1], resol, p1[2], sigma)$fq_new
    scale_new <- 1.2 * scale
    
    # bad fit
    if (fq_new > fit_quality) {
      # optimize scaling factor
      fq <- 0
      scale <- 0
      if (sum(int_vector) > sum(p1[2] * dnorm(mass_vector, p1[1], sigma))) {
        while ((round(fq, digits = 3) != round(fq_new, digits = 3)) && (scale_new < 10000)) {
          fq <- fq_new
          scale <- scale_new
          # fit 1 peak
          fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, use_bounds)
          p1 <- fitted_peak$par
          # get new value for fit quality and scale
          fq_new <- get_fit_quality(mass_vector, int_vector, p1[1], p1[1], resol, p1[2], sigma)$fq_new
          scale_new <- 1.2 * scale
        }
      } else {
        while ((round(fq, digits = 3) != round(fq_new, digits = 3)) && (scale_new < 10000)) {
          fq <- fq_new
          scale <- scale_new
          # fit 1 peak
          fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, use_bounds)
          p1 <- fitted_peak$par
          # get new value for fit quality and scale
          fq_new <- get_fit_quality(mass_vector, int_vector, p1[1], p1[1], resol, p1[2], sigma)$fq_new
          scale_new <- 0.8 * scale
        }
      }
      # use optimized scale factor to fit 1 peak
      if (fq < fq_new) {
        fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, use_bounds)
        p1 <- fitted_peak$par
        fq_new <- fq
      }
    }
  }
  
  roi_value_list <- list("mean" = p1[1], "scale" = p1[2], "sigma" = sigma, "qual" = fq_new)
  return(roi_value_list)
}

fit_2peaks <- function(mass_vector2, mass_vector, int_vector, max_index, resol, use_bounds = FALSE,
                       fit_quality) {
  #' Fit 2 Gaussian peaks in small region of m/z
  #'
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
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
  
  # set range vectors for 2 peaks
  range1 <- c(max_index[1] - 2, max_index[1] - 1, max_index[1], max_index[1] + 1, max_index[1] + 2)
  if (range1[1] == 0) range1 <- range1[-1]
  range2 <- c(max_index[2] - 2, max_index[2] - 1, max_index[2], max_index[2] + 1, max_index[2] + 2)
  if (length(mass_vector) < range2[length(range2)]) range2 <- range2[-length(range2)]
  range1 <- check_overlap(range1, range2)[[1]]
  range2 <- check_overlap(range1, range2)[[2]]
  # check for negative or 0
  remove <- which(range1 < 1)
  if (length(remove) > 0) range1 <- range1[-remove]
  remove <- which(range2 < 1)
  if (length(remove) > 0) range2 <- range2[-remove]
  # remove NA
  if (length(which(is.na(int_vector[range1]))) != 0) range1 <- range1[-which(is.na(int_vector[range1]))]
  if (length(which(is.na(int_vector[range2]))) != 0) range2 <- range2[-which(is.na(int_vector[range2]))]
  
  # fit 2 peaks, first separately, then together
  mu1 <- weighted.mean(mass_vector[range1], int_vector[range1])
  sigma1 <- get_stdev(mass_vector[range1], int_vector[range1])
  fitted_peak <- fit_1gaussian(mass_vector[range1], int_vector[range1], sigma1, mu1, scale, use_bounds)
  p1 <- fitted_peak$par
  # second peak
  mu2 <- weighted.mean(mass_vector[range2], int_vector[range2])
  sigma2 <- get_stdev(mass_vector[range2], int_vector[range2])
  fitted_peak <- fit_1gaussian(mass_vector[range2], int_vector[range2], sigma2, mu2, scale, use_bounds)
  p2 <- fitted_peak$par
  # combined
  fitted_2peaks <- fit_2gaussians(mass_vector, int_vector, sigma1, sigma2, p1[1], p1[2], p2[1], p2[2], use_bounds)
  pc <- fitted_2peaks$par
  
  # get fit quality
  if (is.null(sigma2)) sigma2 <- sigma1
  sum_fit <- (pc[2] * dnorm(mass_vector, pc[1], sigma1)) +
    (pc[4] * dnorm(mass_vector, pc[3], sigma2))
  fq <- get_fit_quality(mass_vector, int_vector, sort(c(pc[1], pc[3]))[1], sort(c(pc[1], pc[3]))[2],
                        resol, sum_fit = sum_fit)$fq_new
  
  # get parameter values
  area1 <- estimate_area(pc[1], resol, pc[2], sigma1, int_factor)
  area2 <- estimate_area(pc[3], resol, pc[4], sigma2, int_factor)
  peak_area <- c(peak_area, area1)
  peak_area <- c(peak_area, area2)
  peak_mean <- c(peak_mean, pc[1])
  peak_mean <- c(peak_mean, pc[3])
  peak_scale <- c(peak_scale, pc[2])
  peak_scale <- c(peak_scale, pc[4])
  peak_sigma <- c(peak_sigma, sigma1)
  peak_sigma <- c(peak_sigma, sigma2)
  
  roi_value_list <- list("mean" = peak_mean, "scale" = peak_scale, "sigma" = peak_sigma, "area" = peak_area, "qual" = fq)
  return(roi_value_list)
}

fit_1gaussian <- function(mass_vector, int_vector, sigma, query_mass, use_bounds) {
  #' Fit a Gaussian curve for a peak with given parameters
  #'
  #' @param mass_vector: Vector of masses (float)
  #' @param int_vector: Vector of intensities (float)
  #' @param sigma: Value for width of the peak (float)
  #' @param query_mass: Value for mass at center of peak (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return opt_fit: list of parameters and values describing the optimal fit
  
  # Initial value used to estimate scaling parameter
  scale <- 2
  
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
    opt_fit <- optim(c(as.numeric(query_mass), as.numeric(scale)),
                     opt_f, control = list(maxit = 10000), method = "L-BFGS-B",
                     lower = lower, upper = upper)
  } else {
    opt_fit <- optim(c(as.numeric(query_mass), as.numeric(scale)),
                     opt_f, control = list(maxit = 10000))
  }
  return(opt_fit)
}

fit_2gaussians <- function(mass_vector, int_vector, sigma1, sigma2,
                           query_mass1, scale1,
                           query_mass2, scale2, use_bounds) {
  #' Fit two Gaussian curves for a peak with given parameters
  #'
  #' @param mass_vector: Vector of masses (float)
  #' @param int_vector: Vector of intensities (float)
  #' @param sigma1: Value for width of the first peak (float)
  #' @param sigma2: Value for width of the second peak (float)
  #' @param query_mass1: Value for mass at center of first peak (float)
  #' @param scale1: Value for scaling intensities for first peak (float)
  #' @param query_mass2: Value for mass at center of second peak (float)
  #' @param scale2: Value for scaling intensities for second peak (float)
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
    if (is.null(query_mass2) && is.null(scale2) && is.null(sigma2)) {
      sigma2 <- sigma1
      opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale1),
                         as.numeric(query_mass1), as.numeric(scale1)),
                       opt_f, control = list(maxit = 10000),
                       method = "L-BFGS-B", lower = lower, upper = upper)
    } else {
      opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale1),
                         as.numeric(query_mass2), as.numeric(scale2)),
                       opt_f, control = list(maxit = 10000),
                       method = "L-BFGS-B", lower = lower, upper = upper)
    }
  } else {
    if (is.null(query_mass2) && is.null(scale2) && is.null(sigma2)) {
      sigma2 <- sigma1
      opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale1),
                         as.numeric(query_mass1), as.numeric(scale1)),
                       opt_f, control = list(maxit = 10000))
    } else {
      opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale1),
                         as.numeric(query_mass2), as.numeric(scale2)),
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

## adapted from getFitQuality.R
# parameter not used: mu_last
get_fit_quality <- function(mass_vector, int_vector, mu_first, mu_last, resol, scale = NULL, sigma = NULL, sum_fit = NULL) {
  #' Fit 1 Gaussian peak in small region of m/z
  #'
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param mu_first: Value for first peak (float)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param sum_fit: Value indicating quality of fit of Gaussian curve (float)
  #'
  #' @return list_params: list of parameters indicating quality of fit (list)
  if (is.null(sum_fit)) {
    mass_vector_int <- mass_vector
    int_vector_int <- int_vector
    # get new fit quality
    fq_new <- mean(abs((scale * dnorm(mass_vector_int, mu_first, sigma)) - int_vector_int) / 
                     rep((max(scale * dnorm(mass_vector_int, mu_first, sigma)) / 2), length(mass_vector_int)))
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

estimate_area <- function(mass_max, resol, scale, sigma) {
  #' Estimate area of Gaussian curve
  #'
  #' @param mass_max: Value for m/z at maximum intensity of a peak (float)
  #' @param resol: Value for resolution (integer)
  #' @param scale: Value for peak width (float)
  #' @param sigma: Value for standard deviation (float)
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #'
  #' @return area_curve: Value for area under the Gaussian curve (float)
  
  # Number used to calculate area under Gaussian curve
  int_factor <- 1 * 10^5 
  
  # avoid vectors that are too big (cannot allocate vector of size ...)
  if (mass_max > 1200) return(0)
  
  # generate a mass_vector with equally spaced m/z values
  fwhm <- get_fwhm(mass_max, resol)
  mz_min <- mass_max - 2 * fwhm
  mz_max <- mass_max + 2 * fwhm
  mz_range <- mz_max - mz_min 
  mass_vector2 <- seq(mz_min, mz_max, length = mz_range * int_factor)
  
  # estimate area under the curve
  area_curve <- sum(scale * dnorm(mass_vector2, mass_max, sigma)) / 100
  
  return(area_curve)  
}

get_fwhm <- function(query_mass, resol) {
  #' Calculate fwhm (full width at half maximum intensity) for a peak
  #'
  #' @param query_mass: Value for mass (float)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return fwhm: Value for full width at half maximum (float)
  
  # set aberrant values of query_mass to zero
  if (is.nan(query_mass)) query_mass <- 0
  if (is.na(query_mass)) query_mass <- 0
  if (is.null(query_mass)) query_mass <- 0
  if (query_mass < 0) query_mass <- 0
  # calculate resolution at given m/z value
  resol_mz <- resol * (1 / sqrt(2) ^ (log2(query_mass / 200)))
  # calculate full width at half maximum
  fwhm <- query_mass / resol_mz
  return(fwhm)
}

fit_optim <- function(mass_vector, int_vector, resol, scanmode) {
  #' Determine optimized fit of Gaussian curve to small region of m/z
  #'
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Vector of intensities for a region of interest (float)
  #' @param resol: Value for resolution (integer)
  #' @param scanmode: Scan mode, positive or negative (string)
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  # Number used to calculate area under Gaussian curve
  int_factor <- 1 * 10^5
  factor <- 1.5
  # Find the index in int_vector with the highest intensity
  max_index <- which(int_vector == max(int_vector))[1]
  mass_max <- mass_vector[max_index]
  int_max <- int_vector[max_index]
  # get peak width
  fwhm <- get_fwhm(mass_max, resol)
  # simplify the peak shape: represent it by a triangle
  mass_max_simple <- c(mass_max - factor * fwhm, mass_max, mass_max + factor * fwhm)
  int_max_simple <- c(0, int_max, 0)
  
  # define mass_diff as difference between last and first value of mass_max_simple
  mass_diff <- mass_max_simple[length(mass_max_simple)] - mass_max_simple[1]
  # generate a second mass_vector with equally spaced m/z values
  mass_vector2 <- seq(mass_max_simple[1], mass_max_simple[length(mass_max_simple)], 
                      length = mass_diff * int_factor)
  sigma <- get_stdev(mass_vector2, int_max_simple)
  # define optimization function for optim based on normal distribution
  opt_f <- function(p, mass_vector, int_vector, sigma, mass_max) {
    curve <- p * dnorm(mass_vector, mass_max, sigma)
    return((max(curve) - max(int_vector))^2)
  }
  
  # get optimal value for fitted Gaussian curve
  opt_fit <- optimize(opt_f, c(0, 100000), tol = 0.0001, mass_vector, int_vector, sigma, mass_max)
  scale <- opt_fit$minimum
  
  # get an estimate of the area under the peak
  area <- estimate_area(mass_max, resol, scale, sigma)
  # put all values for this region of interest into a list
  roi_value_list <- list("mean" = mass_max,
                         "area" = area,
                         "min" = mass_vector2[1],
                         "max" = mass_vector2[length(mass_vector2)])
  return(roi_value_list)
}

within_ppm <- function(mean, scale, sigma, area, mass_vector2, mass_vector, ppm = 4, resol) {
  #' Test whether two mass ranges are within ppm distance of each other
  #'
  #' @param mean: Value for mean m/z (float)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param sigma: Value for standard deviation (float)
  #' @param area: Value for area under the curve (float)
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param ppm: Value for distance between two values of mass (integer)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return list_params: list of parameters indicating quality of fit (list)
  
  # sort
  index <- order(mean)
  mean <- mean[index]
  scale <- scale[index]
  sigma <- sigma[index]
  area <- area[index]
  
  summed <- NULL
  remove <- NULL
  
  if (length(mean) > 1) {
    for (i in 2:length(mean)) {
      if ((abs(mean[i - 1] - mean[i]) / mean[i - 1]) * 10^6 < ppm) {
        
        # avoid double occurrance in sum
        if ((i - 1) %in% summed) next
        
        result_values <- sum_curves(mean[i - 1], mean[i], scale[i - 1], scale[i], sigma[i - 1], sigma[i],
                                    mass_vector2, mass_vector, resol)
        summed <- c(summed, i - 1, i)
        if (is.nan(result_values$mean)) result_values$mean <- 0
        mean[i - 1] <- result_values$mean
        mean[i] <- result_values$mean
        area[i - 1] <- result_values$area
        area[i] <- result_values$area
        scale[i - 1] <- result_values$scale
        scale[i] <- result_values$scale
        sigma[i - 1] <- result_values$sigma
        sigma[i] <- result_values$sigma
        
        remove <- c(remove, i)
      }
    }
  }
  
  if (length(remove) != 0) {
    mean <- mean[-c(remove)]
    area <- area[-c(remove)]
    scale <- scale[-c(remove)]
    sigma <- sigma[-c(remove)]
  }
  
  list_params <- list("mean" = mean, "area" = area, "scale" = scale, "sigma" = sigma, "qual" = NULL)
  return(list_params)
}

sum_curves <- function(mean1, mean2, scale1, scale2, sigma1, sigma2, mass_vector2, mass_vector, resol) {
  #' Sum two curves
  #'
  #' @param mean1: Value for mean m/z of first peak (float)
  #' @param mean2: Value for mean m/z of second peak (float)
  #' @param scale1: Initial value used to estimate scaling parameter for first peak (integer)
  #' @param scale2: Initial value used to estimate scaling parameter for second peak (integer)
  #' @param sigma1: Value for standard deviation for first peak (float)
  #' @param sigma2: Value for standard deviation for second peak (float)
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return list_params: list of parameters indicating quality of fit (list)
  
  sum_fit <- (scale1 * dnorm(mass_vector2, mean1, sigma1)) + (scale2 * dnorm(mass_vector2, mean2, sigma2))
  
  mean1_plus2 <- weighted.mean(c(mean1, mean2), c(max(scale1 * dnorm(mass_vector2, mean1, sigma1)),
                                                  max(scale2 * dnorm(mass_vector2, mean2, sigma2))))
  
  # get new values for parameters
  fwhm <- get_fwhm(mean1_plus2, resol)
  area <- max(sum_fit)
  scale <- scale1 + scale2
  sigma <- (fwhm / 2) * 0.85
  
  list_params <- list("mean" = mean1_plus2, "area" = area, "scale" = scale, "sigma" = sigma)
  return(list_params)
}

