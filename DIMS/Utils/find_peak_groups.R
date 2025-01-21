find_peak_groups <- function(outlist_sorted, mz_tolerance) {
  #' find peaks in all samples with query m/z values and form peak groups
  #'
  #' @param outlist_sorted: matrix of peaks (mzmed, intensity) in all samples
  #' @param mz_tolerance: Value for mass tolerance around query m/z (float)
  #'
  #' @return ints_sorted: matrix of peak groups

  # set up object for intensities for all samples
  ints_allsamps <- matrix(0, nrow = nrow(outlist_sorted), ncol = 3 + (length(sample_names)))
  colnames(ints_allsamps) <- c("mzmed.pgrp", "mzmin.pgrp", "mzmax.pgrp", sample_names)
 
  # start with the m/z with the highest intensity
  row_index <- 1
  while (nrow(outlist_sorted) > 1) {
    # store row numbers
    outlist_sorted$rownr <- 1:nrow(outlist_sorted)
    # find the peaks in the dataset with corresponding m/z plus or minus tolerance
    reference_mass <- outlist_sorted$mzmed.pkt[1]
    minmz_ref <- reference_mass - mz_tolerance
    maxmz_ref <- reference_mass + mz_tolerance
    select_peak_indices <- which((outlist_sorted$mzmed.pkt > minmz_ref) & (outlist_sorted$mzmed.pkt < maxmz_ref))
    select_peaks <- outlist_sorted[select_peak_indices, ]
    # list_of_peaks_used_in_peak_groups_identified <- rbind(list_of_peaks_used_in_peak_groups_identified, tmplist)
    nrsamples <- length(select_peak_indices)
    # if peaks have been found, create a peak group
    if (nrsamples > 0) {
      # if any sample has more than one peak, choose the one with the lowest absolute ppm deviation
      select_peaks$abs_ppmdev <- abs(10^6 * (select_peaks$mzmed.pkt - reference_mass) / reference_mass)
      select_peaks %>% dplyr::group_by(samplenr) %>% dplyr::slice_min(order_by = abs_ppmdev, n = 1)
      # calculate m/z values for peak group
      ints_allsamps[row_index, "mzmed.pgrp"] <- mean(select_peaks$mzmed.pkt)
      ints_allsamps[row_index, "mzmin.pgrp"] <- min(select_peaks$mzmed.pkt)
      ints_allsamps[row_index, "mzmax.pgrp"] <- max(select_peaks$mzmed.pkt)
      # put intensities into proper columns
      column_indices <- c()
      for (sample_name in select_peaks$samplenr) {
        column_number <- which(colnames(ints_allsamps) == sample_name)
        column_indices <- c(column_indices, column_number)
      }
      ints_allsamps[row_index, column_indices] <- select_peaks$height.pkt
      # remove selected peaks from peaklist # Dit gaat niet goed! Indices shiften.
      outlist_sorted <- outlist_sorted[-select_peaks$rownr, ]
      row_index <- row_index + 1
      # print(nrow(outlist_sorted))
    } else {
      outlist_sorted <- outlist_sorted[-1, ]
      row_index <- row_index + 1
    }
  }
 
  # remove empty rows
  ints_allsamps <- ints_allsamps[-which(apply(ints_allsamps, 1, sum) == 0), ]
  # sort by ascending m/z
  ints_allsamps_df <- as.data.frame(ints_allsamps) 
  ints_sorted <- ints_allsamps_df %>% dplyr::arrange(mzmed.pgrp)
 
  return(ints_sorted)
}

