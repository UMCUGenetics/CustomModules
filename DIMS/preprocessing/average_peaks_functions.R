average_peaks_per_sample <- function(peaklist_allrepl_sorted) {
  #' Average the intensity of peaks that occur in different technical replicates of a biological sample
  #'
  #' @param peaklist_allrepl_sorted: Dataframe with peaks sorted on median m/z (float)
  #'
  #' @return averaged_peaks: matrix of averaged peaks (float)
  
  # initialize
  averaged_peaks <- peaklist_allrepl_sorted[0, ]
  
  while (nrow(peaklist_allrepl_sorted) > 1) {
    # store row numbers
    peaklist_allrepl_sorted$rownr <- 1:nrow(peaklist_allrepl_sorted)
    # find the peaks in the dataset with corresponding m/z plus or minus tolerance
    reference_mass <- peaklist_allrepl_sorted$mzmed.pkt[1]
    mz_tolerance <- (reference_mass * ppm_peak) / 10^6
    minmz_ref <- reference_mass - mz_tolerance
    maxmz_ref <- reference_mass + mz_tolerance
    select_peak_indices <- which((peaklist_allrepl_sorted$mzmed.pkt > minmz_ref) & (peaklist_allrepl_sorted$mzmed.pkt < maxmz_ref))
    select_peaks <- peaklist_allrepl_sorted[select_peak_indices, ]
    nrsamples <- length(select_peak_indices)
    # put averaged intensities into a new row
    averaged_1peak <- matrix(0, nrow = 1, ncol = 6) 
    colnames(averaged_1peak) <- c("samplenr", "mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")
    # calculate m/z values for peak group
    averaged_1peak[1, "mzmed.pkt"] <- mean(select_peaks$mzmed.pkt)
    averaged_1peak[1, "mzmin.pkt"] <- min(select_peaks$mzmed.pkt)
    averaged_1peak[1, "mzmax.pkt"] <- max(select_peaks$mzmed.pkt)
    averaged_1peak[1, "fq"] <- nrsamples
    averaged_1peak[1, "height.pkt"] <- mean(select_peaks$height.pkt)
    # remove rownr column and append to averaged_peaks
    peaklist_allrepl_sorted <- peaklist_allrepl_sorted[-select_peaks$rownr, ]
    averaged_peaks <- rbind(averaged_peaks, averaged_1peak)
  }
  # add sample name to first column
  averaged_peaks[ , "samplenr"] <- sample_name
  
  return(averaged_peaks)
}
