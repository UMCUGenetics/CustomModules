# CollectFilled functions

collapse_information <- function(column_label, peakgroup_list, index_dup) {
  #' Collapse identification info for peak groups with the same mass
  #'
  #' @param column_label: Name of column in peakgroup_list (string)
  #' @param peakgroup_list: Peak group list (matrix)
  #' @param index_dup: Index of duplicate peak group (integer)
  #'
  #' @return collapsed_items: Semicolon-separated list of info (string)
  # get the item(s) that need to be collapsed
  list_items <- as.vector(peakgroup_list[index_dup, column_label])
  # remove NA
  if (length(which(is.na(list_items))) > 0) {
    list_items <- list_items[-which(is.na(list_items))]
  }
  collapsed_items <- paste(list_items, collapse = ";")
  return(collapsed_items)
}

merge_duplicate_rows <- function(peakgroup_list) {
  #' Merge identification info for peak groups with the same mass
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #'
  #' @return peakgroup_list_dedup: de-duplicated peak group list (matrix)
  
  options(digits = 16)
  collect <- NULL
  remove <- NULL
  
  # check for peak groups with identical mass
  index_dup <- which(duplicated(peakgroup_list[, "mzmed.pgrp"]))
  
  while (length(index_dup) > 0) {
    # get the index for the peak group which is double
    peaklist_index <- which(peakgroup_list[, "mzmed.pgrp"] == peakgroup_list[index_dup[1], "mzmed.pgrp"])
    single_peakgroup <- peakgroup_list[peaklist_index[1], , drop = FALSE]
    
    # use function collapse_information to concatenate info
    single_peakgroup[, "assi_HMDB"] <- collapse_information("assi_HMDB", peakgroup_list, peaklist_index)
    single_peakgroup[, "iso_HMDB"] <- collapse_information("iso_HMDB", peakgroup_list, peaklist_index)
    single_peakgroup[, "HMDB_code"] <- collapse_information("HMDB_code", peakgroup_list, peaklist_index)
    single_peakgroup[, "all_hmdb_ids"] <- collapse_information("all_hmdb_ids", peakgroup_list, peaklist_index)
    single_peakgroup[, "sec_hmdb_ids"] <- collapse_information("sec_hmdb_ids", peakgroup_list, peaklist_index)
    if (single_peakgroup[, "sec_hmdb_ids"] == ";") single_peakgroup[, "sec_hmdb_ids"] < NA
    
    # keep track of deduplicated entries
    collect <- rbind(collect, single_peakgroup)
    remove <- c(remove, peaklist_index)
    
    # remove current entry from index
    index_dup <- index_dup[-which(peakgroup_list[index_dup, "mzmed.pgrp"] == peakgroup_list[index_dup[1], "mzmed.pgrp"])]
  }
  
  # remove duplicate entries
  if (!is.null(remove)) {
    peakgroup_list <- peakgroup_list[-remove, ]
  }
  # append deduplicated entries
  peakgroup_list_dedup <- rbind(peakgroup_list, collect)
  return(peakgroup_list_dedup)
}

calculate_zscores_peakgrouplist <- function(peakgroup_list) {
  #' Calculate Z-scores for peak groups based on average and standard deviation of controls
  #'
  #' @param peakgroup_list: Peak group list (matrix)
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
  # set all zeros to NA
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
  
  return(peakgroup_list)
}

calculate_ppm_deviation <- function(peakgroup_list) {
  #' Calculate ppm deviation between observed mass and expected theoretical mass
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #'
  #' @return peakgroup_list_ppm: peak group list with ppm column (matrix)
  
  # calculate ppm deviation
  for (row_index in seq_len(nrow(peakgroup_list))) {
    if (!is.na(peakgroup_list$theormz_HMDB[row_index]) &&
        !is.null(peakgroup_list$theormz_HMDB[row_index]) &&
        (peakgroup_list$theormz_HMDB[row_index] != "")) {
      peakgroup_list$ppmdev[row_index] <- 10^6 * (as.numeric(as.vector(peakgroup_list$mzmed.pgrp[row_index])) -
                                                  as.numeric(as.vector(peakgroup_list$theormz_HMDB[row_index]))) /
                                                  as.numeric(as.vector(peakgroup_list$theormz_HMDB[row_index]))
    } else {
      peakgroup_list$ppmdev[row_index] <- NA
    }
  }
  
  return(peakgroup_list)
}

order_columns_peakgrouplist <- function(peakgroup_list) {
  #' Put columns in peak group list in correct order
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #'
  #' @return peakgroup_ordered: peak group list with columns in correct order (matrix)
  
  original_colnames <- colnames(peakgroup_list)
  mass_columns <- c(grep("mzm", original_colnames), grep("nrsamples", original_colnames))
  descriptive_columns <- c(grep("assi_HMDB", original_colnames):grep("avg.int", original_colnames), grep("ppmdev", original_colnames))
  intensity_columns <- c((grep("nrsamples", original_colnames) + 1):(grep("assi_HMDB", original_colnames) - 1))
  # if no Z-scores have been calculated, the following two variables will be empty without consequences for outlist_total
  control_columns <- grep ("ctrls", original_colnames)
  zscore_columns <- grep("_Zscore", original_colnames)
  # create peak group list with columns in correct order
  peakgroup_ordered <- peakgroup_list[ , c(mass_columns, descriptive_columns, intensity_columns, control_columns, zscore_columns)]
  
  return(peakgroup_ordered)
}

