## adapted from sumCurves.R
# variables with fixed values will be removed from function parameters
# plot
# parameter half_max not used
sum_curves <- function(mean1, mean2, scale1, scale2, sigma1, sigma2, mass_vector2, mass_vector, resol, plot) {
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
  #' @param plot: Parameter indicating whether plots should be made (boolean)
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

