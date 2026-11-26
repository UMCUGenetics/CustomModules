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

