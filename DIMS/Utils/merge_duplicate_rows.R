## adapted from mergeDuplicatedRows.R
merge_duplicate_rows <- function(peakgroup_list) {
  #' Merge identification info for peak groups with the same mass
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #'
  #' @return peakgroup_list_dedup: de-duplicated peak group list (matrix)

  collapse <- function(column_label, peakgroup_list, index_dup) {
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
    if (length(which(is.na(list_items))) > 0) list_items <- list_items[-which(is.na(list_items))]
    collapsed_items <- paste(list_items, collapse = ";")
    return(collapsed_items)
  }

  options(digits = 16)
  collect <- NULL
  remove <- NULL

  # check for peak groups with identical mass
  index_dup <- which(duplicated(peakgroup_list[, "mzmed.pgrp"]))

  while (length(index_dup) > 0) {
    # get the index for the peak group which is double
    peaklist_index <- which(peakgroup_list[, "mzmed.pgrp"] == peakgroup_list[index_dup[1], "mzmed.pgrp"])
    single_peakgroup <- peakgroup_list[peaklist_index[1], , drop = FALSE]

    # use function collapse to concatenate info
    single_peakgroup[, "assi_HMDB"] <- collapse("assi_HMDB", peakgroup_list, peaklist_index)
    single_peakgroup[, "iso_HMDB"] <- collapse("iso_HMDB", peakgroup_list, peaklist_index)
    single_peakgroup[, "HMDB_code"] <- collapse("HMDB_code", peakgroup_list, peaklist_index)
    single_peakgroup[, "assi_noise"] <- collapse("assi_noise", peakgroup_list, peaklist_index)
    if (single_peakgroup[, "assi_noise"] == ";") single_peakgroup[, "assi_noise"] <- NA
    single_peakgroup[, "theormz_noise"] <- collapse("theormz_noise", peakgroup_list, peaklist_index)
    if (single_peakgroup[,"theormz_noise"] == "0;0") single_peakgroup[, "theormz_noise"] <- NA

    # keep track of deduplicated entries
    collect <- rbind(collect, single_peakgroup)
    remove <- c(remove, peaklist_index)

    # remove current entry from index
    index_dup <- index_dup[-which(peakgroup_list[index_dup, "mzmed.pgrp"] == peakgroup_list[index_dup[1], "mzmed.pgrp"])]
  }

  # remove duplicate entries
  if (!is.null(remove)) peakgroup_list <- peakgroup_list[-remove, ]
  # append deduplicated entries
  peakgroup_list_dedup <- rbind(peakgroup_list, collect)
  return(peakgroup_list_dedup)
}
