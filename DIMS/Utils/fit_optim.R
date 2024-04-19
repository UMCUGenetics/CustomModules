## adapted from generateGaussian.R
# variables with fixed values will be removed from function parameters
# plot, width, height
# fit_gaussian should be defined before this function is called.
fit_optim <- function(mass_vector, int_vector, resol,
                      plot, scanmode, int_factor, width, height) {
  #' Determine optimized fit of Gaussian curve to small region of m/z
  #'
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Vector of intensities for a region of interest (float)
  #' @param resol: Value for resolution (integer)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param scanmode: Scan mode, positive or negative (string)
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #' @param width: Value for width of plot (integer)
  #' @param height: Value for height of plot (integer)
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)
  
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
                      length = mz_diff * int_factor)
  sigma <- get_stdev(mass_vector2, int_max_simple)
  scale <- optimizeGauss(mass_vector2, int_max_simple, sigma, mass_max)

  # get an estimate of the area under the peak
  area <- getArea(mass_max, resol, scale, sigma, int_factor)
  
  # put all values for this region of interest into a list
  roi_value_list <- list("mean" = mass_max,
                         "area" = area,
                         "min" = mass_vector2[1],
                         "max" = mass_vector2[length(mass_vector2)])
  return(roi_value_list)
}
