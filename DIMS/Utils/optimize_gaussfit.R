## adapted from optimizeGauss.R
optimize_gaussfit <- function(mass_vector, int_vector, sigma, mass_max) {
  #' Optimize fit of Gaussian curve to small region of m/z
  #'
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param int_vector: Vector of intensities for a region of interest (float)
  #' @param sigma: Value for standard deviation (float)
  #' @param mass_max: Value for mass at center of peak (float)
  #'
  #' @return opt_fit: list of fit values for region of interest (list)

  # define optimization function for optim based on normal distribution
  opt_f <- function(p, mass_vector, int_vector, sigma, mass_max) {
    curve <- p * dnorm(mass_vector, mass_max, sigma)
    return((max(curve) - max(int_vector))^2)
  }

  # get optimal value for fitted Gaussian curve
  opt_fit <- optimize(opt_f, c(0, 100000), tol = 0.0001, mass_vector, int_vector, sigma, mass_max)

  return(opt_fit$minimum)
}

