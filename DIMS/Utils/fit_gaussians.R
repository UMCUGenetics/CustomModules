# Gaussian fit functions
## adapted from fitG.R, fit2G.R, fit3G.R and fit4G.R (combined)
fit_1gaussian <- function(mass_vector, int_vector, sigma, query_mass, scale, use_bounds) {
  #' Fit a Gaussian curve for a peak with given parameters
  #'
  #' @param mass_vector: Vector of masses (float)
  #' @param int_vector: Vector of intensities (float)
  #' @param sigma: Value for width of the peak (float)
  #' @param query_mass: Value for mass at center of peak (float)
  #' @param scale: Value for scaling intensities (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return opt_fit: list of parameters and values describing the optimal fit

  # define optimization function for optim based on normal distribution
  opt_f <- function(params) {
    d <- params[2] * dnorm(mass_vector, mean = params[1], sd = sigma)
    sum((d - int_vector) ^ 2)
  }
print("nu in fit_1gaussian")
  if (use_bounds) {
    # determine lower and upper boundaries
    lower <- c(mass_vector[1], 0, mass_vector[1], 0)
    upper <- c(mass_vector[length(mass_vector)], Inf, mass_vector[length(mass_vector)], Inf)
    # get optimal value for fitted Gaussian curve
    tryCatch(opt_fit <- optim(c(as.numeric(query_mass), as.numeric(scale)),
                     opt_f, control = list(maxit = 10000), method = "L-BFGS-B",
                     lower = lower, upper = upper),
	     error = function(e) {
		     # in case of error, use regular optim without boundaries
                     opt_fit <- optim(c(as.numeric(query_mass), as.numeric(scale)),
                                      opt_f, control = list(maxit = 10000))
                     write.table(opt_fit, file = paste0("tryCatch_error_", query_mass, ".txt"), row.names = FALSE)
             } )
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

print("nu in fit_2gaussians")
  if (use_bounds) {
    # determine lower and upper boundaries
    lower <- c(mass_vector[1], 0, mass_vector[1], 0)
    upper <- c(mass_vector[length(mass_vector)], Inf, mass_vector[length(mass_vector)], Inf)
    # get optimal value for 2 fitted Gaussian curves
    if (is.null(query_mass2) && is.null(scale2) && is.null(sigma2)) {
      sigma2 <- sigma1
      tryCatch(opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale1),
                                as.numeric(query_mass1), as.numeric(scale1)),
                                opt_f, control = list(maxit = 10000),
                                method = "L-BFGS-B", lower = lower, upper = upper),
	       error = function(e) {
                     # in case of error, use regular optim without boundaries
                     opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale1),
                                      as.numeric(query_mass1), as.numeric(scale1)),
                                      opt_f, control = list(maxit = 10000))
                     write.table(res, file = paste0("tryCatch_error_2gauss_", query_mass, ".txt"), row.names = FALSE) } )
    } else {
      tryCatch(opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale1),
                                as.numeric(query_mass2), as.numeric(scale2)),
                                opt_f, control = list(maxit = 10000),
                                method = "L-BFGS-B", lower = lower, upper = upper),
	       error = function(e) {
                     # in case of error, use regular optim without boundaries
                     opt_fit <- optim(c(as.numeric(query_mass1), as.numeric(scale1),
                                      as.numeric(query_mass1), as.numeric(scale1)),
                                      opt_f, control = list(maxit = 10000))
                     write.table(res, file = paste0("tryCatch_error_2gauss_else_", query_mass, ".txt"), row.names = FALSE) } )
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


fit_3gaussians <- function(mass_vector, int_vector, sigma1, sigma2, sigma3,
                           query_mass1, scale1,
                           query_mass2, scale2,
                           query_mass3, scale3, use_bounds) {
  #' Fit three Gaussian curves for a peak with given parameters
  #'
  #' @param mass_vector: Vector of masses (float)
  #' @param int_vector: Vector of intensities (float)
  #' @param sigma1: Value for width of the first peak (float)
  #' @param sigma2: Value for width of the second peak (float)
  #' @param sigma3: Value for width of the third peak (float)
  #' @param query_mass1: Value for mass at center of first peak (float)
  #' @param scale1: Value for scaling intensities for first peak (float)
  #' @param query_mass2: Value for mass at center of second peak (float)
  #' @param scale2: Value for scaling intensities for second peak (float)
  #' @param query_mass3: Value for mass at center of third peak (float)
  #' @param scale3: Value for scaling intensities for third peak (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return opt_fit: list of parameters and values describing the optimal fit

  # define optimization function for optim based on normal distribution
  opt_f <- function(params) {
    d <- params[2] * dnorm(mass_vector, mean = params[1], sd = sigma1) +
         params[4] * dnorm(mass_vector, mean = params[3], sd = sigma2) +
         params[6] * dnorm(mass_vector, mean = params[5], sd = sigma3)
    sum((d - int_vector) ^ 2)
  }

print("nu in fit_3gaussians")
  if (use_bounds) {
    # determine lower and upper boundaries
    lower <- c(mass_vector[1], 0, mass_vector[1], 0, mass_vector[1], 0)
    upper <- c(mass_vector[length(mass_vector)], Inf, mass_vector[length(mass_vector)], Inf,
               mass_vector[length(mass_vector)], Inf)
    # get optimal value for 3 fitted Gaussian curves
    opt_fit <- optim(c(query_mass1, scale1,
                       query_mass2, scale2,
                       query_mass3, scale3),
                     opt_f, control = list(maxit = 10000),
                     method = "L-BFGS-B", lower = lower, upper = upper)
  } else {
    opt_fit <- optim(c(query_mass1, scale1,
                       query_mass2, scale2,
                       query_mass3, scale3),
                     opt_f, control = list(maxit = 10000))
  }
  return(opt_fit)
}

fit_4gaussians <- function(mass_vector, int_vector, sigma1, sigma2, sigma3, sigma4,
                           query_mass1, scale1,
                           query_mass2, scale2,
                           query_mass3, scale3,
                           query_mass4, scale4, use_bounds) {
  #' Fit four Gaussian curves for a peak with given parameters
  #'
  #' @param mass_vector: Vector of masses (float)
  #' @param int_vector: Vector of intensities (float)
  #' @param sigma1: Value for width of the first peak (float)
  #' @param sigma2: Value for width of the second peak (float)
  #' @param sigma3: Value for width of the third peak (float)
  #' @param sigma4: Value for width of the fourth peak (float)
  #' @param query_mass1: Value for mass at center of first peak (float)
  #' @param scale1: Value for scaling intensities for first peak (float)
  #' @param query_mass2: Value for mass at center of second peak (float)
  #' @param scale2: Value for scaling intensities for second peak (float)
  #' @param query_mass3: Value for mass at center of third peak (float)
  #' @param scale3: Value for scaling intensities for third peak (float)
  #' @param query_mass4: Value for mass at center of fourth peak (float)
  #' @param scale4: Value for scaling intensities for fourth peak (float)
  #' @param use_bounds: Boolean to indicate whether boundaries are to be used
  #'
  #' @return opt_fit: list of parameters and values describing the optimal fit

  # define optimization function for optim based on normal distribution
  opt_f <- function(params) {
    d <- params[2] * dnorm(mass_vector, mean = params[1], sd = sigma1) +
         params[4] * dnorm(mass_vector, mean = params[3], sd = sigma2) +
         params[6] * dnorm(mass_vector, mean = params[5], sd = sigma3) +
         params[8] * dnorm(mass_vector, mean = params[7], sd = sigma4)
    sum((d - int_vector) ^ 2)
  }

print("nu in fit_4gaussians")
  if (use_bounds) {
    # determine lower and upper boundaries
    lower <- c(mass_vector[1], 0, mass_vector[1], 0, mass_vector[1], 0, mass_vector[1], 0)
    upper <- c(mass_vector[length(mass_vector)], Inf, mass_vector[length(mass_vector)], Inf,
               mass_vector[length(mass_vector)], Inf, mass_vector[length(mass_vector)], Inf)
    # get optimal value for 3 fitted Gaussian curves
    opt_fit <- optim(c(query_mass1, scale1,
                       query_mass2, scale2,
                       query_mass3, scale3,
                       query_mass4, scale4),
                     opt_f, control = list(maxit = 10000),
                     method = "L-BFGS-B", lower = lower, upper = upper)
  } else {
    opt_fit <- optim(c(query_mass1, scale1,
                       query_mass2, scale2,
                       query_mass3, scale3,
                       query_mass4, scale4),
                     opt_f, control = list(maxit = 10000))
  }
  return(opt_fit)
}
