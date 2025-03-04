## adapted from fit1Peak.R, fit2peaks.R, fit3peaks.R and fit4peaks.R (combined)
# variables with fixed values will be removed from function parameters
# plot, int_factor
fit_1peak <- function(mass_vector2, mass_vector, int_vector, max_index, scale, resol, plot, fit_quality, use_bounds) {
  #' Fit 1 Gaussian peak in small region of m/z
  #'
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param max_index: Index in int_vector with the highest intensity (integer)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param fit_quality: Value indicating quality of fit of Gaussian curve (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  if (length(int_vector) < 3) {
    message("Range too small, no fit possible!")
  } else {
    if ((length(int_vector) == 4)) {
      # fit 1 peak
      mu <- weighted.mean(mass_vector, int_vector)
      sigma <- get_stdev(mass_vector, int_vector)
      fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, scale, use_bounds)
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
      #### bug fix ####
      sigma <- get_stdev(mass_vector, int_vector)
      if (sum(int_vector[range1]) == 0) return(list("mean" = mean(mass_vector), "scale" = 1, "sigma" = sigma, "qual" = 0))
 
      # remove NA
      if (length(which(is.na(int_vector[range1]))) != 0) {
        range1 <- range1[-which(is.na(int_vector[range1]))]
      }
      # fit 1 peak
      mu <- weighted.mean(mass_vector[range1], int_vector[range1])
      sigma <- get_stdev(mass_vector[range1], int_vector[range1])
      fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, scale, use_bounds)
    }
    
    p1 <- fitted_peak$par
    print(paste0("mass_vector:", mass_vector, "int_vector:", int_vector, "sigma:", sigma, "p1[1]:", p1[1], "p1[2]:", p1[2]))
    if (is.null(p1)) {
     roi_value_list <- list("mean" = 0, "scale" = 0, "sigma" = 0, "area" = 0, "qual" = 0)
     return(roi_value_list)
    }
  
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
          fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, scale, use_bounds)
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
          fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, scale, use_bounds)
          p1 <- fitted_peak$par
          # get new value for fit quality and scale
          fq_new <- get_fit_quality(mass_vector, int_vector, p1[1], p1[1], resol, p1[2], sigma)$fq_new
          scale_new <- 0.8 * scale
        }
      }
      # use optimized scale factor to fit 1 peak
      if (fq < fq_new) {
        fitted_peak <- fit_1gaussian(mass_vector, int_vector, sigma, mu, scale, use_bounds)
        p1 <- fitted_peak$par
        fq_new <- fq
      }
    }
  }
 
  roi_value_list <- list("mean" = p1[1], "scale" = p1[2], "sigma" = sigma, "qual" = fq_new)
  return(roi_value_list)	
}

fit_2peaks <- function(mass_vector2, mass_vector, int_vector, max_index, scale, resol, use_bounds = FALSE,
                       plot = FALSE, fit_quality, int_factor) {
  #' Fit 2 Gaussian peaks in small region of m/z
  #'
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param max_index: Index in int_vector with the highest intensity (integer)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param fit_quality: Value indicating quality of fit of Gaussian curve (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
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
  if (is.null(pc)) {
    roi_value_list <- list("mean" = 0, "scale" = 0, "sigma" = 0, "area" = 0, "qual" = 0)
    return(roi_value_list)
  }
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

fit_3peaks <- function(mass_vector2, mass_vector, int_vector, max_index, scale, resol, use_bounds = FALSE, 
                       plot = FALSE, fit_quality, int_factor) {
  #' Fit 3 Gaussian peaks in small region of m/z
  #'
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param max_index: Index in int_vector with the highest intensity (integer)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param fit_quality: Value indicating quality of fit of Gaussian curve (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  peak_mean <- NULL
  peak_area <- NULL
  peak_scale <- NULL
  peak_sigma <- NULL
  
  # set range vectors for 3 peaks
  range1 <- c(max_index[1] - 2, max_index[1] - 1, max_index[1], max_index[1] + 1, max_index[1] + 2)
  range2 <- c(max_index[2] - 2, max_index[2] - 1, max_index[2], max_index[2] + 1, max_index[2] + 2)
  range3 <- c(max_index[3] - 2, max_index[3] - 1, max_index[3], max_index[3] + 1, max_index[3] + 2)
  remove <- which(range1 < 1)
  if (length(remove) > 0) {
    range1 <- range1[-remove]
  }
  remove <- which(range2 < 1)
  if (length(remove) > 0) {
    range2 <- range2[-remove]
  }
  if (length(mass_vector) < range3[length(range3)]) range3 <- range3[-length(range3)]
  range1 <- check_overlap(range1, range2)[[1]]
  range2 <- check_overlap(range1, range2)[[2]]
  range2 <- check_overlap(range2, range3)[[1]]
  range3 <- check_overlap(range2, range3)[[2]]
  # check for negative or 0
  remove <- which(range1 < 1)
  if (length(remove) > 0) range1 <- range1[-remove]
  remove <- which(range2 < 1)
  if (length(remove) > 0) range2 <- range2[-remove]
  remove <- which(range3 < 1)
  if (length(remove) > 0) range3 <- range3[-remove]
  # remove NA
  if (length(which(is.na(int_vector[range1]))) != 0) range1 <- range1[-which(is.na(int_vector[range1]))]
  if (length(which(is.na(int_vector[range2]))) != 0) range2 <- range2[-which(is.na(int_vector[range2]))]
  if (length(which(is.na(int_vector[range3]))) != 0) range3 <- range3[-which(is.na(int_vector[range3]))]
  
  # fit 3 peaks, first separately, then together
  mu1 <- weighted.mean(mass_vector[range1], int_vector[range1])
  sigma1 <- get_stdev(mass_vector[range1], int_vector[range1])
  fitted_peak <- fit_1gaussian(mass_vector[range1], int_vector[range1], sigma1, mu1, scale, use_bounds)
  p1 <- fitted_peak$par
  # second peak
  mu2 <- weighted.mean(mass_vector[range2], int_vector[range2])
  sigma2 <- get_stdev(mass_vector[range2], int_vector[range2])
  fitted_peak <- fit_1gaussian(mass_vector[range2], int_vector[range2], sigma2, mu2, scale, use_bounds)
  p2 <- fitted_peak$par
  # third peak
  mu3 <- weighted.mean(mass_vector[range3], int_vector[range3])
  sigma3 <- get_stdev(mass_vector[range3], int_vector[range3])
  fitted_peak <- fit_1gaussian(mass_vector[range3], int_vector[range3], sigma3, mu3, scale, use_bounds)
  p3 <- fitted_peak$par
  # combined
  fitted_3peaks <- fit_3gaussians(mass_vector, int_vector, sigma1, sigma2, sigma3,
                                  p1[1], p1[2], p2[1], p2[2], p3[1], p3[2], use_bounds)
  pc <- fitted_3peaks$par
  print(pc)
  if (is.null(pc)) {
    roi_value_list <- list("mean" = 0, "scale" = 0, "sigma" = 0, "area" = 0, "qual" = 0)
    return(roi_value_list)
  }
  # get fit quality
  sum_fit = (pc[2] * dnorm(mass_vector, pc[1], sigma1)) +
    (pc[4] * dnorm(mass_vector, pc[3], sigma2)) +
    (pc[6] * dnorm(mass_vector, pc[5], sigma3))
  fq <- get_fit_quality(mass_vector, int_vector, sort(c(pc[1], pc[3], pc[5]))[1], sort(c(pc[1], pc[3], pc[5]))[3],
                        resol, sum_fit = sum_fit)$fq_new
  
  # get parameter values
  area1 <- estimate_area(pc[1], resol, pc[2], sigma1, int_factor)
  area2 <- estimate_area(pc[3], resol, pc[4], sigma2, int_factor)
  area3 <- estimate_area(pc[5], resol, pc[6], sigma3, int_factor)
  peak_area <- c(peak_area, area1)
  peak_area <- c(peak_area, area2)
  peak_area <- c(peak_area, area3)
  peak_mean <- c(peak_mean, pc[1])
  peak_mean <- c(peak_mean, pc[3])
  peak_mean <- c(peak_mean, pc[5])
  peak_scale <- c(peak_scale, pc[2])
  peak_scale <- c(peak_scale, pc[4])
  peak_scale <- c(peak_scale, pc[6])
  peak_sigma <- c(peak_sigma, sigma1)
  peak_sigma <- c(peak_sigma, sigma2)
  peak_sigma <- c(peak_sigma, sigma3)
  
  roi_value_list <- list("mean" = peak_mean, "scale" = peak_scale, "sigma" = peak_sigma, "area" = peak_area, "qual" = fq)
  return(roi_value_list)
}

fit_4peaks <- function(mass_vector2, mass_vector, int_vector, max_index, scale, resol, use_bounds = FALSE,
                       plot = FALSE, fit_quality, int_factor) {
  #' Fit 4 Gaussian peaks in small region of m/z
  #'
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param max_index: Index in int_vector with the highest intensity (integer)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param fit_quality: Value indicating quality of fit of Gaussian curve (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
  peak_mean <- NULL
  peak_area <- NULL
  peak_scale <- NULL
  peak_sigma <- NULL
  
  # set range vectors for 4 peaks
  range1 <- c(max_index[1] - 2, max_index[1] - 1, max_index[1], max_index[1] + 1, max_index[1] + 2)
  range2 <- c(max_index[2] - 2, max_index[2] - 1, max_index[2], max_index[2] + 1, max_index[2] + 2)
  range3 <- c(max_index[3] - 2, max_index[3] - 1, max_index[3], max_index[3] + 1, max_index[3] + 2)
  range4 <- c(max_index[4] - 2, max_index[4] - 1, max_index[4], max_index[4] + 1, max_index[4] + 2)
  if (range1[1] == 0) range1 <- range1[-1]
  if (length(mass_vector) < range4[length(range4)]) range4 <- range4[-length(range4)]
  range1 <- check_overlap(range1, range2)[[1]]
  range2 <- check_overlap(range1, range2)[[2]]
  range2 <- check_overlap(range2, range3)[[1]]
  range3 <- check_overlap(range2, range3)[[2]]
  range3 <- check_overlap(range3, range4)[[1]]
  range4 <- check_overlap(range3, range4)[[2]]
  remove <- which(range4 > length(mass_vector))
  if (length(remove) > 0) {
    range4 <- range4[-remove]
  }
  # check for negative or 0
  remove <- which(range1 < 1)
  if (length(remove) > 0) range1 <- range1[-remove]
  remove <- which(range2 < 1)
  if (length(remove) > 0) range2 <- range2[-remove]
  remove <- which(range3 < 1)
  if (length(remove) > 0) range3 <- range3[-remove]
  remove <- which(range4 < 1)
  if (length(remove) > 0) range4 <- range4[-remove]
  # remove NA
  if (length(which(is.na(int_vector[range1]))) != 0) range1 <- range1[-which(is.na(int_vector[range1]))]
  if (length(which(is.na(int_vector[range2]))) != 0) range2 <- range2[-which(is.na(int_vector[range2]))]
  if (length(which(is.na(int_vector[range3]))) != 0) range3 <- range3[-which(is.na(int_vector[range3]))]
  if (length(which(is.na(int_vector[range4]))) != 0) range4 <- range4[-which(is.na(int_vector[range4]))]
  
  # fit 4 peaks, first separately, then together
  mu1 <- weighted.mean(mass_vector[range1], int_vector[range1])
  sigma1 <- get_stdev(mass_vector[range1], int_vector[range1])
  fitted_peak <- fit_1gaussian(mass_vector[range1], int_vector[range1], sigma1, mu1, scale, use_bounds)
  p1 <- fitted_peak$par
  # second peak
  mu2 <- weighted.mean(mass_vector[range2], int_vector[range2])
  sigma2 <- get_stdev(mass_vector[range2], int_vector[range2])
  fitted_peak <- fit_1gaussian(mass_vector[range2], int_vector[range2], sigma2, mu2, scale, use_bounds)
  p2 <- fitted_peak$par
  # third peak
  mu3 <- weighted.mean(mass_vector[range3], int_vector[range3])
  sigma3 <- get_stdev(mass_vector[range3], int_vector[range3])
  fitted_peak <- fit_1gaussian(mass_vector[range3], int_vector[range3], sigma3, mu3, scale, use_bounds)
  p3 <- fitted_peak$par
  # fourth peak
  mu4 <- weighted.mean(mass_vector[range4], int_vector[range4])
  sigma4 <- get_stdev(mass_vector[range4], int_vector[range4])
  fitted_peak <- fit_1gaussian(mass_vector[range4], int_vector[range4], sigma4, mu4, scale, use_bounds)
  p4 <- fitted_peak$par
  # combined
  fitted_4peaks <- fit_4gaussians(mass_vector, int_vector, sigma1, sigma2, sigma3, sigma3, 
                                  p1[1], p1[2], p2[1], p2[2], p3[1], p3[2],  p4[1], p4[2], use_bounds)
  pc <- fitted_4peaks$par
  print(pc)
  if (is.null(pc)) {
    roi_value_list <- list("mean" = 0, "scale" = 0, "sigma" = 0, "area" = 0, "qual" = 0)
    return(roi_value_list)
  }
  
  # get fit quality
  sum_fit <-  (pc[2] * dnorm(mass_vector, pc[1], sigma1)) + 
    (pc[4] * dnorm(mass_vector, pc[3], sigma2)) + 
    (pc[6] * dnorm(mass_vector, pc[5], sigma3)) + 
    (pc[8] * dnorm(mass_vector, pc[7], sigma3))
  fq <- get_fit_quality(mass_vector, int_vector, 
                        sort(c(pc[1], pc[3], pc[5], pc[7]))[1], sort(c(pc[1], pc[3], pc[5], pc[7]))[4], 
                        resol, sum_fit = sum_fit)$fq_new
  
  # get parameter values
  area1 <- estimate_area(pc[1], resol, pc[2], sigma1, int_factor)
  area2 <- estimate_area(pc[3], resol, pc[4], sigma2, int_factor)
  area3 <- estimate_area(pc[5], resol, pc[6], sigma3, int_factor)
  area4 <- estimate_area(pc[7], resol, pc[8], sigma4, int_factor)
  peak_area <- c(peak_area, area1)
  peak_area <- c(peak_area, area2)
  peak_area <- c(peak_area, area3)
  peak_area <- c(peak_area, area4)
  peak_mean <- c(peak_mean, pc[1])
  peak_mean <- c(peak_mean, pc[3])
  peak_mean <- c(peak_mean, pc[5])
  peak_mean <- c(peak_mean, pc[7])
  peak_scale <- c(peak_scale, pc[2])
  peak_scale <- c(peak_scale, pc[4])
  peak_scale <- c(peak_scale, pc[6])
  peak_scale <- c(peak_scale, pc[8])
  peak_sigma <- c(peak_sigma, sigma1)
  peak_sigma <- c(peak_sigma, sigma2)
  peak_sigma <- c(peak_sigma, sigma3)
  peak_sigma <- c(peak_sigma, sigma4)
  
  roi_value_list <- list("mean" = peak_mean, "scale" = peak_scale, "sigma" = peak_sigma, "area" = peak_area, "qual" = fq)
  return(roi_value_list)
}

