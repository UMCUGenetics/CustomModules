sum_intensities_adducts <- function(peakgroup_list, hmdb_part, adducts, z_score) {
  #' Sum intensities for different adducts of the same metabolite
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #' @param hmdb_part: Matrix of metabolites , part of the HMDB (matrix)
  #' @param adducts: Vector of adducts (vector of integers)
  #' @param z_score: Value indicating whether Z-scores have been calculated (integer)
  #'
  #' @return adductsum: peak group list with summed intensities (matrix)
  hmdb_part_info <- cbind(HMDB_id = rownames(hmdb_part), CompoundName = hmdb_part[, "CompoundName"])

  # create overview of row indices for each metabolite_adduct combination in peaklist
  # split the all_hmdb_ids column into list with each id as a value
  hmdb_in_peaklist <- strsplit(peakgroup_list$all_hmdb_ids, ";")
  # avoid rows with only "" in HMDB_code column
  hmdb_in_peaklist[which(hmdb_in_peaklist == "")] <- ";"
  hmdb_in_peaklist_rownr <- c()

  # create dataframe with for each HMDB id a row number
  hmdb_in_peaklist_rownr <- data.frame(
    row_id = rep(seq_along(hmdb_in_peaklist), lengths(hmdb_in_peaklist)),
    hmdb_id = unlist(hmdb_in_peaklist)
  )
  # remove empty rows and duplicates
  hmdb_in_peaklist_rownr <- hmdb_in_peaklist_rownr %>%
    filter(hmdb_id != "") %>%
    distinct()

  # find intensity columns in peakgroup_list
  if (z_score == 1) {
    int_cols_ctrls <- grep("C", colnames(peakgroup_list)[1:which(colnames(peakgroup_list) == "avg.ctrls")])
    int_cols_pats <- grep("P", colnames(peakgroup_list)[1:which(colnames(peakgroup_list) == "avg.ctrls")])
    int_cols <- c(int_cols_ctrls, int_cols_pats)
  } else {
    int_cols_start <- which(colnames(peakgroup_list) == "nrsamples") + 1
    int_cols_end <- which(colnames(peakgroup_list) == "assi_HMDB") - 1
    int_cols <- c(int_cols_start:int_cols_end)
  }

  # initialize
  names <- NULL
  adductsum <- NULL
  names_long <- NULL

  # find adducts of each metabolite and sum the intensities
  if (nrow(hmdb_part_info) == 0) {
    return(adductsum)
  }

  for (hmdb_index in seq_len(nrow(hmdb_part_info))) {
    compound <- hmdb_part_info[hmdb_index, "HMDB_id"]
    compound_plus_adducts <- c(compound, paste(compound, adducts, sep = "_"))

    # find indices of rows in peakgroup_list that contain compound plus adducts
    metab_row <- which(hmdb_in_peaklist_rownr$hmdb_id %in% compound_plus_adducts)
    metab_indices <- as.numeric(hmdb_in_peaklist_rownr$row_id[metab_row])

    # find intensities and sum them
    ints <- peakgroup_list[metab_indices, int_cols]
    total <- apply(ints, 2, sum)

    # add to adductsum
    if (sum(total) != 0) {
      names <- c(names, compound)
      adductsum <- rbind(adductsum, total)
      names_long <- c(names_long, hmdb_part_info[hmdb_index, "CompoundName"])
    }
  }

  if (!is.null(adductsum)) {
    rownames(adductsum) <- names
    adductsum <- cbind(adductsum, "HMDB_name" = names_long)
    # Add HMDB info
    cols_hmdb_info <- c("HMDB_ID_all", "sec_HMDB_ID", "HMDB_name_all, theormz_HMDB")
    hmdb_info <- hmdb_part[names, cols_hmdb_info, drop = FALSE]
    adductsum <- cbind(adductsum, hmdb_info)
  }

  return(adductsum)
}
