## adapted from fitGaussianInit.R
# variables with fixed values will be removed from function parameters
# int_factor, scale, outdir, plot, thresh, width, height
# mz_index, start_index, end_index, sample_name not used.
# fit_gaussian should be defined before this function is called.
fit_init <- function(mass_vector, int_vector, int_factor, scale, resol,
                     outdir, sample_name, scanmode, plot, width, height,
                     mz_index, start_index, end_index) {
  #' Initial fit of Gaussian curve to small region of m/z
  #'
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Vector of intensities for a region of interest (float)
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param outdir: Path for output directory (string)
  #' @param sample_name: Sample name (string)
  #' @param scanmode: Scan mode, positive or negative (string)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param width: Value for width of plot (integer)
  #' @param height: Value for height of plot (integer)
  #' @param mz_index: Index of m/z value with non-zero intensity (integer)
  #' @param start_index: Index of start of m/z range in mass_vector (integer)
  #' @param end_index: Index of end of m/z range in mass_vector (integer)
  #'
  #' @return roi_value_list: list of fit values for region of interest (list)

  # define mass_diff as difference between last and first value of mass_vector
  mass_diff <- mass_vector[length(mass_vector)] - mass_vector[1]
  # generate a second mass_vector with equally spaced m/z values
  mass_vector2 <- seq(mass_vector[1], mass_vector[length(mass_vector)],
                      length = mass_diff * int_factor)

  # Find the index in int_vector with the highest intensity
  max_index <- which(int_vector == max(int_vector))

  roi_values <- fit_gaussian(mass_vector2, mass_vector, int_vector, max_index,
                             scale, resol, outdir, force = length(max_index),
                             useBounds = FALSE, plot, scanmode, int_factor, width, height)

  roi_value_list <- list("mean" = roi_values$mean,
                         "area" = roi_values$area,
                         "qual" = roi_values$qual,
                         "min" = roi_values$min,
                         "max" = roi_values$max)

  return(roi_value_list)
}

