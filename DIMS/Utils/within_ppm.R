## adapted from isWithinXppm.R
# variables with fixed values will be removed from function parameters
# plot
within_ppm <- function(mean, scale, sigma, area, mass_vector2, mass_vector, ppm = 4, resol, plot) {
  #' Test whether two mass ranges are within ppm distance of each other
  #'
  #' @param mean: Value for mean m/z (float)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param sigma: Value for standard deviation (float)
  #' @param area: Value for area under the curve (float)
  #' @param mass_vector2: Vector of equally spaced m/z values (float)
  #' @param mass_vector: Vector of m/z values for a region of interest (float)
  #' @param ppm: Value for distance between two values of mass (integer)
  #' @param resol: Value for resolution (integer)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #'
  #' @return list_params: list of parameters indicating quality of fit (list)

  # sort
  index <- order(mean)
  mean <- mean[index]
  scale <- scale[index]
  sigma <- sigma[index]
  area <- area[index]

  summed <- NULL
  remove <- NULL

  if (length(mean) > 1) {
    for (i in 2:length(mean)) {
      if ((abs(mean[i - 1] - mean[i]) / mean[i - 1]) * 10^6 < ppm) {

        # avoid double occurance in sum
        if ((i - 1) %in% summed) next

        result_values <- sum_curves(mean[i - 1], mean[i], scale[i - 1], scale[i], sigma[i - 1], sigma[i],
                                    mass_vector2, mass_vector, resol, plot)
        summed <- c(summed, i - 1, i)
        if (is.nan(result_values$mean)) result_values$mean <- 0
        mean[i - 1] <- result_values$mean
        mean[i] <- result_values$mean
        area[i - 1] <- result_values$area
        area[i] <- result_values$area
        scale[i - 1] <- result_values$scale
        scale[i] <- result_values$scale
        sigma[i - 1] <- result_values$sigma
        sigma[i] <- result_values$sigma

        remove <- c(remove, i)
      }
    }
  }

  if (length(remove) != 0) {
    mean <- mean[-c(remove)]
    area <- area[-c(remove)]
    scale <- scale[-c(remove)]
    sigma <- sigma[-c(remove)]
  }

  list_params <- list("mean" = mean, "area" = area, "scale" = scale, "sigma" = sigma, "qual" = NULL)
  return(list_params)
}

