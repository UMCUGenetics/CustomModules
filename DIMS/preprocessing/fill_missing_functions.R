fill_missing_intensities <- function(peakgroup_list, repl_pattern, thresh) {
  #' Replace intensities that are zero with random value
  #'
  #' @param peakgroup_list: Peak groups (matrix)
  #' @param repl_pattern: Replication pattern (list of strings)
  #' @param thresh: Value for threshold between noise and signal (integer)
  #'
  #' @return final_outlist: peak groups with filled-in intensities (matrix)

  # replace missing intensities with random values around threshold
  if (!is.null(peakgroup_list)) {
    for (sample_index in 1:length(names(repl_pattern))) {
      sample_peaks <- peakgroup_list[, names(repl_pattern)[sample_index]]
      zero_intensity <- which(sample_peaks <= 0)
      if (!length(zero_intensity)) {
        next
      }
      for (zero_index in 1:length(zero_intensity)) {
        peakgroup_list[zero_intensity[zero_index], names(repl_pattern)[sample_index]] <- rnorm(n = 1, 
                                                                                               mean = thresh,
                                                                                               sd = 100)
      }
    }

    # Add column with average intensity; find intensity columns first
    int_cols <- which(colnames(peakgroup_list) %in% names(repl_pattern))
    peakgroup_list <- cbind(peakgroup_list, "avg.int" = apply(peakgroup_list[, int_cols], 1, mean))

    return(peakgroup_list)
  }
}
