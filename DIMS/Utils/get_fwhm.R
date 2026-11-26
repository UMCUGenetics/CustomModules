## adapted from getFwhm.R
get_fwhm <- function(query_mass, resol) {
  #' Calculate fwhm (full width at half maximum intensity) for a peak
  #'
  #' @param query_mass: Value for mass (float)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return fwhm: Value for full width at half maximum (float)

  # set aberrant values of query_mass to zero
  if (is.nan(query_mass)) query_mass <- 0
  if (is.na(query_mass)) query_mass <- 0
  if (is.null(query_mass)) query_mass <- 0
  if (query_mass < 0) query_mass <- 0
  # calculate resolution at given m/z value
  resol_mz <- resol * (1 / sqrt(2) ^ (log2(query_mass / 200)))
  # calculate full width at half maximum
  fwhm <- query_mass / resol_mz
  return(fwhm)
}

