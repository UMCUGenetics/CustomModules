estimate_area <- function(mass_max, resol, scale, sigma, int_factor) {
  #' Estimate area of Gaussian curve
  #'
  #' @param mass_max: Value for m/z at maximum intensity of a peak (float)
  #' @param resol: Value for resolution (integer)
  #' @param scale: Value for peak width (float)
  #' @param sigma: Value for standard deviation (float)
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #'
  #' @return area_curve: Value for area under the Gaussian curve (float)

  # generate a mass_vector with equally spaced m/z values
  fwhm <- get_fwhm(mass_max, resol)
  mz_min <- mass_max - 2 * fwhm
  mz_max <- mass_max + 2 * fwhm
  mz_range <- mz_max - mz_min 
  mass_vector2 <- seq(mz_min, mz_max, length = 1000)

  # estimate area under the curve
  area_curve <- sum(scale * dnorm(mass_vector2, mass_max, sigma)) / 100

  return(area_curve)  
}

