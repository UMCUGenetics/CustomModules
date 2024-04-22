## adapted from searchMZRange.R
# variables with fixed values will be removed from function parameters
# int_factor, scale, outdir, plot, thresh, width, height
# allpeaks_values should be generated here, not passed on from do_peakfinding
search_mzrange <- function(ints_fullrange, allpeaks_values, int_factor, scale, resol,
                           outdir, sample_name, scanmode, plot, width, height, thresh) {
  #' Divide the full m/z range into regions of interest with min, max and mean m/z
  #'
  #' @param ints_fullrange: Named list of intensities (float)
  #' @param allpeaks_values: Empty list to store results for all peaks
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param outdir: Path for output directory (string)
  #' @param sample_name: Sample name (string)
  #' @param scanmode: Scan mode, positive or negative (string)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param width: Value for width of plot (integer)
  #' @param height: Value for height of plot (integer)
  #' @param thresh: Value for noise level threshold (integer)
  #'
  #' @return allpeaks_values: list of m/z regions of interest

  # find indices where intensity is not equal to zero
  nonzero_indices <- as.vector(which(ints_fullrange != 0))

  # bad infusion. These should have been taken out in AverageTechReplicates
  if (length(nonzero_indices) == 0) return(allpeaks_values)

  # initialize
  end_index <- NULL
  start_index <- nonzero_indices[1]
  # maximum length of region of interest
  max_roi_length <- 15

  # find regions of interest
  for (mz_index in 1:length(nonzero_indices)) {
    # check whether mz_index is smaller than length(nonzero_indices).
    # only false if mz_index == length(nonzero_indixes). 
    # second check is true at the end of a peak.
    if (mz_index < length(nonzero_indices) && (nonzero_indices[mz_index + 1] - nonzero_indices[mz_index]) > 1) {
      end_index <- nonzero_indices[mz_index]
      # get m/z values and intensities for this region of interest
      mass_vector <- as.numeric(names(ints_fullrange)[c(start_index:end_index)])
      int_vector <- as.vector(ints_fullrange[c(start_index:end_index)])
      # check whether the vector of intensities is not empty. 
      if (length(int_vector) != 0) {
        # check if intensity is above threshold or the maximum intensity is NaN
        if (max(int_vector) < thresh || is.nan(max(int_vector))) {
          # go to next region of interest
          start_index <- nonzero_indices[mz_index + 1]
          next
        }
        # check if there are more intensities than maximum for region of interest
        if (length(int_vector) > max_roi_length) {
          # trim lowest intensities to zero
          int_vector[which(int_vector < min(int_vector) * 1.1)] <- 0
          # split the range into multiple sub ranges
          sub_range <- int_vector
          names(sub_range) <- mass_vector
          allpeaks_values <- search_mzrange(sub_range, allpeaks_values, int_factor,
                                            scale, resol, outdir, sample_name, scanmode,
                                            plot, width, height, thresh)
        # A proper peak needs to have at least 3 intensities above threshold  
        } else if (length(int_vector) > 3) {
          # check if the sum of intensities is above zero. Why is this necessary?
          if (sum(int_vector) == 0) next
          # get initial fit values
          roi_values <- fit_init(mass_vector, int_vector, int_factor, scale, resol,
                                 outdir, sample_name, scanmode, plot, width, height,
                                 mz_index, start_index, end_index)

          if (roi_values$qual[1] == 1) {
            # get optimized fit values
            roi_values <- fit_optim(mass_vector, int_vector, resol, plot,
                                    scanmode, int_factor, width, height)
            # add region of interest to list of all peaks
            allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean)
            allpeaks_values$area <- c(allpeaks_values$area, roi_values$area)
            allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
            allpeaks_values$min <- c(allpeaks_values$min, roi_values$min)
            allpeaks_values$max <- c(allpeaks_values$max, roi_values$max)
            allpeaks_values$qual <- c(allpeaks_values$qual, 0)
            allpeaks_values$spikes <- allpeaks_values$spikes + 1

          } else {
            for (j in 1:length(roi_values$mean)){
              allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean[j])
              allpeaks_values$area <- c(allpeaks_values$area, roi_values$area[j])
              allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
              allpeaks_values$min <- c(allpeaks_values$min, roi_values$min[1])
              allpeaks_values$max <- c(allpeaks_values$max, roi_values$max[1])
              allpeaks_values$qual <- c(allpeaks_values$qual, roi_values$qual[1])
            }
          }

        } else {

          roi_values <- fit_optim(mass_vector, int_vector, resol,
                                  plot, scanmode, int_factor, width, height)

          allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean)
          allpeaks_values$area <- c(allpeaks_values$area, roi_values$area)
          allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
          allpeaks_values$min <- c(allpeaks_values$min, roi_values$min)
          allpeaks_values$max <- c(allpeaks_values$max, roi_values$max)
          allpeaks_values$qual <- c(allpeaks_values$qual, 0)
          allpeaks_values$spikes <- allpeaks_values$spikes + 1
        }
      }
      start_index <- nonzero_indices[mz_index + 1]
    }
  }

  # last little range
  end_index <- nonzero_indices[length(nonzero_indices)]
  mass_vector <- as.numeric(names(ints_fullrange)[c(start_index:end_index)])
  int_vector <- as.vector(ints_fullrange[c(start_index:end_index)])

  if (length(int_vector) != 0) {
    # check if intensity above threshold
    if (max(int_vector) < thresh || is.nan(max(int_vector))) {
      # do nothing
    } else {
      # check if there are more intensities than maximum for region of interest
      if (length(int_vector) > max_roi_length) {
        # trim lowest intensities to zero
        int_vector[which(int_vector < min(int_vector) * 1.1)] <- 0
        # split the range into multiple sub ranges
        sub_range <- int_vector
        names(sub_range) <- mass_vector

        allpeaks_values <- search_mzrange(sub_range, allpeaks_values, int_factor, scale, resol,
                                          outdir, sample_name, scanmode,
                                          plot, width, height, thresh)

      } else if (length(int_vector) > 3) {
        # Check only zeros
        if (sum(int_vector) == 0) next

        roi_values <- fit_init(mass_vector, int_vector, int_factor, scale,
                               resol, outdir, sample_name, scanmode,
                               plot, width, height,
                               mz_index, start_index, end_index)

        if (roi_values$qual[1] == 1) {
          roi_values <- fit_optim(mass_vector, int_vector, resol,
                                  plot, scanmode, int_factor, width, height)

          allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean)
          allpeaks_values$area <- c(allpeaks_values$area, roi_values$area)
          allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
          allpeaks_values$min <- c(allpeaks_values$min, roi_values$min)
          allpeaks_values$max <- c(allpeaks_values$max, roi_values$max)
          allpeaks_values$qual <- c(allpeaks_values$qual, 0)
          allpeaks_values$spikes <- allpeaks_values$spikes + 1

        } else {
          for (j in 1:length(roi_values$mean)){
            allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean[j])
            allpeaks_values$area <- c(allpeaks_values$area, roi_values$area[j])
            allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
            allpeaks_values$min <- c(allpeaks_values$min, roi_values$min[1])
            allpeaks_values$max <- c(allpeaks_values$max, roi_values$max[1])
            allpeaks_values$qual <- c(allpeaks_values$qual, roi_values$qual[1])
          }
        }
      } else {
        roi_values <- fit_optim(mass_vector, int_vector, resol,
                                plot, scanmode, int_factor, width, height)

        allpeaks_values$mean <- c(allpeaks_values$mean, roi_values$mean)
        allpeaks_values$area <- c(allpeaks_values$area, roi_values$area)
        allpeaks_values$nr <- c(allpeaks_values$nr, sample_name)
        allpeaks_values$min <- c(allpeaks_values$min, roi_values$min)
        allpeaks_values$max <- c(allpeaks_values$max, roi_values$max)
        allpeaks_values$qual <- c(allpeaks_values$qual, 0)
        allpeaks_values$spikes <- allpeaks_values$spikes + 1
      }
    }
  }
  return(allpeaks_values)
}

