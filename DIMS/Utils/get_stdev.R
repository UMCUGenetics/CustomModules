## adapted from getSD.R
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

