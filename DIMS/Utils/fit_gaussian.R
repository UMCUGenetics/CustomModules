## adapted from fitGaussian.R
# variables with fixed values will be removed from function parameters
# scale, outdir, plot, width, height
# max_index doesn't need to be passed to this function, can be determined here.
# remove plot sections (commented out)
# several functions need to be loaded before this function can run
fit_gaussian <- function(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                        outdir, force, use_bounds, plot, scanmode,
                        int_factor, width, height) {
  #' Fit 1, 2, 3 or 4 Gaussian peaks in small region of m/z
  #'
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Value used to calculate area under Gaussian curve (integer)
  #' @param max_index: Index in int_vector with the highest intensity (integer)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param outdir: Path for output directory (string)
  #' @param force: Number of local maxima in int_vector (integer)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param scanmode: Scan mode, positive or negative (string)
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #' @param width: Value for width of plot (integer)
  #' @param height: Value for height of plot (integer)
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)

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
    fit_values <- fit_1peak(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                            plot, fit_quality1, use_bounds)
    # set initial value for scale factor
    scale <- 2
    # test if the mean is outside the m/z range
    if (fit_values$mean[1] < mass_vector[1] || fit_values$mean[1] > mass_vector[length(mass_vector)]) {
      # run this function again with fixed boundaries
      return(fit_gaussian(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                          outdir, force = 1, use_bounds = TRUE, plot, scanmode, int_factor, width, height))

    } else {
      # test if the fit is bad
      if (fit_values$qual > fit_quality1) {
        # Try to fit two curves; find two local maxima
        new_index <- which(diff(sign(diff(int_vector))) == -2) + 1
        # test if there are two indices in new_index
        if (length(new_index) != 2) {
          new_index <- round(length(mass_vector) / 3)
          new_index <- c(new_index, 2 * new_index)
        }
        # run this function again with two local maxima
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, new_index,
                            scale, resol, outdir, force = 2, use_bounds = FALSE,
                            plot, scanmode, int_factor, width, height))
      # good fit
      } else {
        peak_mean <- c(peak_mean, fit_values$mean)
        peak_area <- c(peak_area, getArea(fit_values$mean, resol, fit_values$scale,
                                          fit_values$sigma, int_factor))
        peak_qual <- fit_values$qual
        peak_min <- mass_vector[1]
        peak_max <- mass_vector[length(mass_vector)]
      }
    }

  #### Two local maxima; need at least 6 data points for this ####
  } else if (force == 2 && (length(mass_vector) > 6)) {
    # determine fit values for 2 Gaussian peaks (mean, scale, sigma, qual)
    fit_values <- fit_2peaks(mass_vector2, mass_vector, int_vector, new_index, scale, resol,
                             use_bounds, plot, fit_quality, int_factor)
    # test if one of the means is outside the m/z range
    if (fit_values$mean[1] < mass_vector[1] || fit_values$mean[1] > mass_vector[length(mass_vector)] ||
        fit_values$mean[2] < mass_vector[1] || fit_values$mean[2] > mass_vector[length(mass_vector)]) {
      # check if fit quality is bad
      if (fit_values$qual > fit_quality) {
        # run this function again with fixed boundaries
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                            outdir, force = 2, use_bounds = TRUE,
                            plot, scanmode, int_factor, width, height))
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
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, new_index,
                            scale, resol, outdir, force = 3, use_bounds = FALSE,
                            plot, scanmode, int_factor, width, height))
      # good fit, all means are within m/z range
      } else {
        # check if means are within 3 ppm and sum if so
        tmp <- fit_values$qual
        nr_means_new <- -1
        nr_means <- length(fit_values$mean)
        while (nr_means != nr_means_new) {
          nr_means <- length(fit_values$mean)
          fit_values <- within_ppm(fit_values$mean, fit_values$scale, fit_values$sigma, fit_values$area, 
                                     mass_vector2, mass_vector, ppm = 4, resol, plot)
          nr_means_new <- length(fit_values$mean)
        }
        # restore original quality score
        fit_values$qual <- tmp
      }  
    }

  # Three local maxima; need at least 6 data points for this 
  } else if (force == 3 && (length(mass_vector) > 6)) {
    # determine fit values for 3 Gaussian peaks (mean, scale, sigma, qual)
    fit_values <- fit_3peaks(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                             use_bounds, plot, fit_quality, int_factor)
    # test if one of the means is outside the m/z range
    if (fit_values$mean[1] < mass_vector[1] || fit_values$mean[1] > mass_vector[length(mass_vector)] ||
        fit_values$mean[2] < mass_vector[1] || fit_values$mean[2] > mass_vector[length(mass_vector)] ||
        fit_values$mean[3] < mass_vector[1] || fit_values$mean[3] > mass_vector[length(mass_vector)]) {

      # check if fit quality is bad
      if (fit_values$qual > fit_quality) {
        # run this function again with fixed boundaries
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                            outdir, force, use_bounds = TRUE,
                            plot, scanmode, int_factor, width, height))
      } else {
        # check which mean is outside range and remove it from the list of means
        # NB: peak_mean and other variables have not been given values from 2-peak fit yet!
        for (i in 1:length(fit_values$mean)) {
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
        # Try to fit four curves; find four local maxima
        new_index <- which(diff(sign(diff(int_vector))) == -2) + 1
        # test if there are four indices in new_index
        if (length(new_index) != 4) {
          new_index <- round(length(mass_vector) / 5)
          new_index <- c(new_index, 2 * new_index, 3 * new_index, 4 * new_index)
        }
        # run this function again with four local maxima
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, new_index, scale, resol,
                            outdir, force = 4, use_bounds = FALSE, plot, scanmode,
                            int_factor, width, height))
      # good fit, all means are within m/z range
      } else {
        # check if means are within 4 ppm and sum if so  
        tmp <- fit_values$qual
        nr_means_new <- -1
        nr_means <- length(fit_values$mean)
        while (nr_means != nr_means_new) {
          nr_means <- length(fit_values$mean)
          fit_values <- within_ppm(fit_values$mean, fit_values$scale, fit_values$sigma, fit_values$area, 
                                   mass_vector2, mass_vector, ppm = 4, resol, plot)
          nr_means_new <- length(fit_values$mean)
        }
        # restore original quality score
        fit_values$qual <- tmp
      }
    }
    
  #### Four local maxima; need at least 6 data points for this ####  
  } else if (force == 4 && (length(mass_vector) > 6)) {
    # determine fit values for 4 Gaussian peaks (mean, scale, sigma, qual)
    fit_values <- fit_4peaks(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                             use_bounds, plot, fit_quality, int_factor)
    # test if one of the means is outside the m/z range
    if (fit_values$mean[1] < mass_vector[1] || fit_values$mean[1] > mass_vector[length(mass_vector)] ||
        fit_values$mean[2] < mass_vector[1] || fit_values$mean[2] > mass_vector[length(mass_vector)] ||
        fit_values$mean[3] < mass_vector[1] || fit_values$mean[3] > mass_vector[length(mass_vector)] ||
        fit_values$mean[4] < mass_vector[1] || fit_values$mean[4] > mass_vector[length(mass_vector)]) {

      # check if quality of fit is bad
      if (fit_values$qual > fit_quality) {
        # run this function again with fixed boundaries
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                            outdir, force, use_bounds = TRUE,
                            plot, scanmode, int_factor, width, height))

      } else {
        # check which mean is outside range and remove it from the list of means
        # NB: peak_mean and other variables have not been given values from 2-peak fit yet!
        for (i in 1:length(fit_values$mean)) {
          if (fit_values$mean[i] < mass_vector[1] | fit_values$mean[i] > mass_vector[length(mass_vector)]) {
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
        # Try to fit 1 curve, force = 5
        return(fit_gaussian(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                            outdir, force = 5, use_bounds = FALSE,
                            plot, scanmode, int_factor, width, height))
      # good fit, all means are within m/z range
      } else {
        # check if means are within 4 ppm and sum if so
        tmp <- fit_values$qual 
        nr_means_new <- -1
        nr_means <- length(fit_values$mean)
        while (nr_means != nr_means_new) {
          nr_means <- length(fit_values$mean)
          fit_values <- within_ppm(fit_values$mean, fit_values$scale, fit_values$sigma, fit_values$area, 
                                   mass_vector2, mass_vector, ppm = 4, resol, plot)
          nr_means_new <- length(fit_values$mean)
        }
        # restore original quality score
        fit_values$qual <- tmp
      }  
    }
    
  # More than four local maxima; fit 1 peak.
  } else {
    scale <- 2
    fit_quality1 <- 0.40
    use_bounds <- TRUE
    max_index <- which(int_vector == max(int_vector))
    fit_values <- fit_1peak(mass_vector2, mass_vector, int_vector, max_index, scale, resol,
                            plot, fit_quality1, use_bounds)
    # check for bad fit
    if (fit_values$qual > fit_quality1) {
      # remove
      if (plot) dev.off()    
      # get fit values from fit_optim
      fit_values <- fit_optim(mass_vector, int_vector, resol, plot, scanmode, int_factor, width, height)
      peak_mean <- c(peak_mean, fit_values$mean)
      peak_area <- c(peak_area, fit_values$area)
      peak_min <- fit_values$min
      peak_max <- fit_values$max
      peak_qual <- 0
    } else {
      peak_mean <- c(peak_mean, fit_values$mean)
      peak_area <- c(peak_area, get_area(fit_values$mean, resol, fit_values$scale, fit_values$sigma, int_factor))
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

