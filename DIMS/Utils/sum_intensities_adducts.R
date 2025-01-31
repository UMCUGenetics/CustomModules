sum_intensities_adducts <- function(peakgroup_list, hmdb_part, adducts, z_score) {
  #' Sum intensities for different adducts of the same metabolite
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #' @param hmdb_part: Matrix of metabolites , part of the HMDB (matrix)
  #' @param adducts: Vector of adducts (vector of integers)
  #' @param z_score: Value indicating whether Z-scores have been calculated (integer)
  #'
  #' @return adductsum: peak group list with summed intensities (matrix)
  hmdb_codes <- rownames(hmdb_part)
  hmdb_names <- hmdb_part[, 1]

  # create overview of row indices for each metabolite_adduct combination in peaklist
  hmdb_in_peaklist <- peakgroup_list$HMDB_code
  # avoid rows with only "" in HMDB_code column
  hmdb_in_peaklist[which(hmdb_in_peaklist == "")] <- ";"
  hmdb_in_peaklist_rownr <- c()
  for (row_nr in 1:length(hmdb_in_peaklist)) {
    hmdb_split <- strsplit(hmdb_in_peaklist[row_nr], ";")[[1]]
    hmdb_rownr <- cbind(hmdb_split, row_nr)
    hmdb_in_peaklist_rownr <- rbind(hmdb_in_peaklist_rownr, hmdb_rownr)
  }
  hmdb_in_peaklist_rownr <- as.data.frame(hmdb_in_peaklist_rownr)
  # remove NA, if any
  if (sum(is.na(hmdb_in_peaklist_rownr$hmdb_split)) > 0 ) {
    hmdb_in_peaklist_rownr <- hmdb_in_peaklist_rownr[-which(is.na(hmdb_in_peaklist_rownr$hmdb_split)), ]
  }
  
  # find intensity columns in peakgroup_list
  if (z_score == 1) {
    int_cols_C <- grep("C", colnames(peakgroup_list)[1:which(colnames(peakgroup_list) == "avg.ctrls")])
    int_cols_P <- grep("P", colnames(peakgroup_list)[1:which(colnames(peakgroup_list) == "avg.ctrls")])
    int_cols <- c(int_cols_C, int_cols_P)
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
  if (length(hmdb_codes) == 0) {
    return(adductsum)
  }

  for (hmdb_index in 1:length(hmdb_codes)) {
    compound <- hmdb_codes[hmdb_index]
    compound_plus <- c(compound, paste(compound, adducts, sep = "_"))

    # find indices of rows in peakgroup_list that contain compound plus adducts
    metab_row <- which(hmdb_in_peaklist_rownr$hmdb_split %in% compound_plus)
    metab_indices <- as.numeric(hmdb_in_peaklist_rownr$row_nr[metab_row])
      
    # find intensities and sum them
    ints <- peakgroup_list[metab_indices, int_cols]
    total <- apply(ints, 2, sum)
      
    # add to adductsum
    if (sum(total) != 0) {
      names <- c(names, compound)
      adductsum <- rbind(adductsum, total)
      names_long <- c(names_long, hmdb_names[hmdb_index])
    }
  }
    
  if (!is.null(adductsum)) {
    rownames(adductsum) <- names
    adductsum <- cbind(adductsum, "HMDB_name" = names_long)
  } 
 
  return(adductsum)
}
