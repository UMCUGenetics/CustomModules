# functions for peak finding

search_regions_of_interest <- function(ints_fullrange) {
  #' Divide the full m/z range into regions of interest with indices
  #'
  #' @param ints_fullrange: Matrix with m/z values and intensities (float)
  #'
  #' @return regions_of_interest: matrix of m/z regions of interest (integer)

  # find indices where intensity is not equal to zero
  nonzero_positions <- as.vector(which(ints_fullrange$int != 0))

  # find regions of interest (look for consecutive numbers)
  regions_of_interest_consec <- seqToIntervals(nonzero_positions)
  # add length of regions of interest
  regions_of_interest_diff <- regions_of_interest_consec[, 2] - regions_of_interest_consec[, 1] + 1
  regions_of_interest_length <- cbind(regions_of_interest_consec, length = regions_of_interest_diff)
  # remove short lengths; a peak should have at least 3 data points
  if (any(regions_of_interest_length[, "length"] < 3)) {
    regions_of_interest_gte3 <- regions_of_interest_length[-which(regions_of_interest_length[, "length"] < 3), ]
  } else {
    regions_of_interest_gte3 <- regions_of_interest_length
  }
  # test for length of roi (region of interest). If length is greater than 11, break up into separate rois
  remove_roi_index <- c()
  new_rois_all <- regions_of_interest_gte3[0, ]
  for (roi_nr in 1:nrow(regions_of_interest_gte3)) {
    if (regions_of_interest_gte3[roi_nr, "length"] > 11) {
      roi <- ints_fullrange[(regions_of_interest_gte3[roi_nr, "from"]:regions_of_interest_gte3[roi_nr, "to"]), ]
      roi_intrange <- as.numeric(roi$int)
      roi_firstindex <- as.numeric(rownames(roi)[1])
      # look for local minima that separate the peaks
      local_min_positions <- which(diff(sign(diff(roi_intrange))) == 2) + 1
      if (length(local_min_positions) > 0) {
        remove_roi_index <- c(remove_roi_index, roi_nr)
        # find new indices for rois after splitting
        new_rois_splitroi <- as.data.frame(matrix(0, ncol = 3, nrow = (length(local_min_positions) + 1)))
        colnames(new_rois_splitroi) <- colnames(regions_of_interest_gte3)
	# fill new rois matrix; from in column 1, to in column 2 and length in column 3
        new_rois_splitroi[, 1] <- c(roi_firstindex, roi_firstindex + local_min_positions)
        new_rois_splitroi[, 2] <- c(roi_firstindex + local_min_positions, roi_firstindex + length(roi_intrange))
        new_rois_splitroi[, 3] <- new_rois_splitroi[, 2] - new_rois_splitroi[, 1]
        # append
        new_rois_all <- rbind(new_rois_all, new_rois_splitroi)
      } else {
        # if there are no local minima, all intensities belong to one peak.
        remove_roi_index <- c(remove_roi_index, roi_nr)
        # look for maximum and take 5 intensities to the left and right
        max_index <- which(roi_intrange == max(roi_intrange))
        max_pos <- as.numeric(rownames(roi)[max_index])
        new_roi <- data.frame(from = 0, to = 0, length = 0)
        new_roi[, 1] <- max_pos - 5
        new_roi[, 2] <- max_pos + 5
        new_roi[, 3] <- new_roi[, 2] - new_roi[, 1] + 1
        # append
        new_rois_all <- rbind(new_rois_all, new_roi)
      }
    }
  }

  # remove rois that have been split into chunks or shortened
  if (length(remove_roi_index) > 0) {
    regions_of_interest_minus_short <- regions_of_interest_gte3[-remove_roi_index, ]
  } else {
    regions_of_interest_minus_short <- regions_of_interest_gte3
  }
  # combine remaining rois with info on chunks
  regions_of_interest_split <- rbind(regions_of_interest_minus_short, new_rois_all)
  # remove regions of interest with short lengths again
  if (any(regions_of_interest_split[, "length"] < 3)) {
    regions_of_interest_final <- regions_of_interest_split[-which(regions_of_interest_split[, "length"] < 3), ]
  } else {
    regions_of_interest_final <- regions_of_interest_split
  }

  # sort on first index
  if (nrow(regions_of_interest_final) > 1){
    regions_of_interest_sorted <- regions_of_interest_final %>% as.data.frame %>% dplyr::arrange(from)
  } else {
    regions_of_interest_sorted <- regions_of_interest_final
  }

  return(regions_of_interest_sorted)
}

integrate_peaks <- function(ints_fullrange, regions_of_interest, resol, peak_thresh) {
  #' Fit Gaussian peak for each region of interest and integrate area under the curve
  #'
  #' @param ints_fullrange: Named list of intensities (float)
  #' @param regions_of_interest: Named list of intensities (float)
  #' @param resol: Value for resolution (integer)
  #' @param peak_thresh: Value for noise level threshold (integer) # NOT USED YET
  #'
  #' @return allpeaks_values: matrix of integrated peaks

  # initialize dataframe to store results for all peaks
  allpeaks_values <- matrix(0, nrow = nrow(regions_of_interest), ncol = 5)
  colnames(allpeaks_values) <- c("mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")

  for (roi_nr in 1:nrow(regions_of_interest)) {
    # find m/z values and intensities corresponding to the region of interest
    index_range <- regions_of_interest[roi_nr, "from"] : regions_of_interest[roi_nr, "to"]
    roi_mzrange <- as.numeric(ints_fullrange$mz[index_range])
    roi_intrange <- as.numeric(ints_fullrange$int[index_range])
    # find m/z value for maximum of peak (mu in Gaussian function)
    weighted_mzmean <- weighted.mean(roi_mzrange, roi_intrange)
    # find expected peak width at this m/z value
    fwhm_mzmed <- get_fwhm(weighted_mzmean, resol)
    # calculate sigma for Gaussian curve (https://en.wikipedia.org/wiki/Full_width_at_half_maximum)
    sigma_peak <- fwhm_mzmed / 2.355
    # find scale factor for Gaussian. Initial estimate is maximum intensity in roi
    # divided by 2 for better correlation with intensities from old PeakFinding method
    scale_factor <- max(roi_intrange) / 2
    # fit Gaussian peak. Find intensities according to normal distribution
    normal_ints <- scale_factor * gaussfunc(roi_mzrange, weighted_mzmean, sigma_peak)
    # sum intensities in roi
    sum_ints <- sum(roi_intrange)
    sum_gauss <- sum(normal_ints)
    # calculate quality of fit
    quality_fit <- 1 - sum(abs(normal_ints - roi_intrange)) / sum_ints
    # put all values into dataframe
    allpeaks_values[roi_nr, "mzmed.pkt"] <- weighted_mzmean
    allpeaks_values[roi_nr, "fq"] <- quality_fit
    allpeaks_values[roi_nr, "mzmax.pkt"] <- max(roi_mzrange)
    allpeaks_values[roi_nr, "mzmin.pkt"] <- min(roi_mzrange)
    allpeaks_values[roi_nr, "height.pkt"] <- sum_gauss
  }

  # remove peaks with height = 0, look for NA or NaN
  remove_na <- which(is.na(as.numeric(allpeaks_values[, "height.pkt"])))
  if (length(remove_na) > 0) {
    allpeaks_values <- allpeaks_values[-remove_na, ]
  }

  return(allpeaks_values)
}

# in the next version of the docker image, package R.utils will be included
seqToIntervals <- function(idx) {
  #' Find consecutive stretches of numbers
  #' function seqToIntervals copied from R.utils library
  #' see https://rdrr.io/cran/R.utils/src/R/seqToIntervals.R
  #'
  #' @param idx: Sequence of indices (integers)
  #'
  #' @return res: Matrix of start and end positions of consecutive numbers (matrix)

  # Clean up sequence
  idx <- as.integer(idx)
  idx <- unique(idx)
  idx <- sort(idx)

  len_idx <- length(idx)
  if (len_idx == 0L) {
    res <- matrix(NA_integer_, nrow = 0L, ncol = 2L)
    colnames(res) <- c("from", "to")
    return(res)
  }

  # Identify end points of intervals
  diff_idx <- diff(idx)
  diff_idx <- (diff_idx > 1)
  diff_idx <- which(diff_idx)
  nr_intervals <- length(diff_idx) + 1

  # Allocate return matrix
  res <- matrix(NA_integer_, nrow = nr_intervals, ncol = 2L)
  colnames(res) <- c("from", "to")

  from_value <- idx[1]
  to_value <- from_value - 1
  last_value <- from_value

  count <- 1
  for (running_index in seq_along(idx)) {
    value <- idx[running_index]
    if (value - last_value > 1) {
      to_value <- last_value
      res[count, ] <- c(from_value, to_value)
      from_value <- value
      count <- count + 1
    }
    last_value <- value
  }

  if (to_value < from_value) {
    to_value <- last_value
    res[count, ] <- c(from_value, to_value)
  }

  return(res)
}

get_fwhm <- function(query_mass, resol) {
  #' Calculate fwhm (full width at half maximum intensity) for a peak
  #'
  #' @param query_mass: Value for mass (float)
  #' @param resol: Value for resolution (integer)
  #'
  #' @return fwhm: Value for full width at half maximum (float)

  # set aberrant values of query_mass to default value of 200
  if (is.na(query_mass)) {
    query_mass <- 200
  }
  if (query_mass < 0) {
    query_mass <- 200
  }
  # calculate resolution at given m/z value
  resol_mz <- resol * (1 / sqrt(2) ^ (log2(query_mass / 200)))
  # calculate full width at half maximum
  fwhm <- query_mass / resol_mz

  return(fwhm)
}


# from https://rdrr.io/cran/rvmethod/src/R/gaussfit.R
#' Gaussian Function from package rvmethod
#'
#' This function returns the unnormalized (height of 1.0) Gaussian curve with a
#' given center and spread.
#'
#' @param x the vector of values at which to evaluate the Gaussian
#' @param mu the center of the Gaussian
#' @param sigma the spread of the Gaussian (must be greater than 0)
#' @return vector of values of the Gaussian
#' @examples x = seq(-4, 4, length.out = 100)
#' y = gaussfunc(x, 0, 1)
#' plot(x, y)
#'
#' @import stats
#'
#' @export
gaussfunc <- function(x, mu, sigma) {
  return(exp(-((x - mu) ^ 2) / (2 * (sigma ^ 2))))
}

