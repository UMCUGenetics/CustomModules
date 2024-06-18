## adapted from statistics_z.R
# refactor: change column names from avg.ctrls to avg_ctrls, sd.ctrls to sd_ctrls
# check logic of parameter adducts
calculate_zscores <- function(peakgroup_list, adducts) {
  #' Calculate Z-scores for peak groups based on average and standard deviation of controls
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #' @param sort_col: Column to sort on (string)
  #' @param adducts: Parameter indicating whether there are adducts in the list (boolean)
  #'
  #' @return peakgroup_list_dedup: de-duplicated peak group list (matrix)

  case_label <- "P"
  control_label <- "C"
  # get index for new column names
  startcol <- ncol(peakgroup_list) + 3

  # calculate mean and standard deviation for Control group
  ctrl_cols <- grep(control_label, colnames(peakgroup_list), fixed = TRUE)
  case_cols <- grep(case_label, colnames(peakgroup_list), fixed = TRUE)
  int_cols <- c(ctrl_cols, case_cols)
  # set al zeros to NA
  peakgroup_list[, int_cols][peakgroup_list[, int_cols] == 0] <- NA
  ctrl_ints <- peakgroup_list[, ctrl_cols, drop = FALSE]
  peakgroup_list$avg.ctrls <- apply(ctrl_ints, 1, function(x) mean(as.numeric(x), na.rm = TRUE))
  peakgroup_list$sd.ctrls <- apply(ctrl_ints, 1, function(x) sd(as.numeric(x), na.rm = TRUE))

  # set new column names and calculate Z-scores
  colnames_zscores <- NULL
  for (col_index in int_cols) {
    col_name <- colnames(peakgroup_list)[col_index]
    colnames_zscores <- c(colnames_zscores, paste0(col_name, "_Zscore"))
    zscores_1col <- (as.numeric(as.vector(unlist(peakgroup_list[, col_index]))) -
                     peakgroup_list$avg.ctrls) / peakgroup_list$sd.ctrls
    peakgroup_list <- cbind(peakgroup_list, zscores_1col)
  }

  # apply new column names to columns at end plus avg and sd columns
  colnames(peakgroup_list)[startcol:ncol(peakgroup_list)] <- colnames_zscores

  # add ppm deviation column
  zscore_cols <- grep("Zscore", colnames(peakgroup_list), fixed = TRUE)
  if (!adducts) {
    if ((dim(peakgroup_list[, zscore_cols])[2] + 6) != (startcol - 1)) {
      ppmdev <- array(1:nrow(peakgroup_list), dim = c(nrow(peakgroup_list)))
      # calculate ppm deviation
      for (i in 1:nrow(peakgroup_list)) {
        if (!is.na(peakgroup_list$theormz_HMDB[i]) && 
            !is.null(peakgroup_list$theormz_HMDB[i]) &&
            (peakgroup_list$theormz_HMDB[i] != "")) {
          ppmdev[i] <- 10^6 * (as.numeric(as.vector(peakgroup_list$mzmed.pgrp[i])) -
                               as.numeric(as.vector(peakgroup_list$theormz_HMDB[i]))) /
                               as.numeric(as.vector(peakgroup_list$theormz_HMDB[i]))
        } else {
          ppmdev[i] <- NA
        }
      }
      peakgroup_list <- cbind(peakgroup_list[, 1:6], ppmdev = ppmdev, peakgroup_list[, 7:ncol(peakgroup_list)])
    }
  }

  return(peakgroup_list)
}
