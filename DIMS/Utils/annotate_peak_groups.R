annotate_peak_groups <- function(ints_sorted, hmdb_add_iso) {
  #' annotate peak groups; assign  metabolites (adducts, isotopes) with suitable mass from HMDB
  #'
  #' @param ints_sorted: matrix of peak groups
  #' @param hmdb_add_iso: subset of HMDB (matrix)
  #'
  #' @return peakgrouplist_identified: matrix of peak groups with annotation

  # Initialize matrix for annotation
  assigned_hmdb <- matrix("", nrow = nrow(ints_sorted), ncol = 7)
  colnames(assigned_hmdb) <- c("assi_HMDB", "all_hmdb_names", "iso_HMDB", "HMDB_code",
                               "all_hmdb_ids", "sec_hmdb_ids", "theormz_HMDB")
 
  # for each peak group, find all entries in HMDB part with mass within ppm range
  for (row_number in 1:nrow(ints_sorted)) {
    # initialize to make sure there's no information from the previous peak group
    select_metabolites <- NULL
    select_adducts <- NULL
    select_isotopes <- NULL
    all_hmdb_names <- ""
    all_hmdb_ids <- ""
    all_isotope_names <- ""
    first_hmdb_name <- ""
    first_hmdb_code <- ""
    sec_hmdb_ids <- ""
    theor_mz <- 0
 
    # take reference mass
    reference_mass <- ints_sorted[row_number, "mzmed.pgrp"]
    # select indices for all HMDB entries with mass between +/- ppm tolerance
    select_from_hmdb <- which(hmdb_add_iso[ , column_label] > (reference_mass - mz_tolerance) & 
                                hmdb_add_iso[ , column_label] < (reference_mass + mz_tolerance))
    if (length(select_from_hmdb) > 0) {
      # get dataframe of all entries which are selected
      select_hmdb_df <- hmdb_add_iso[select_from_hmdb, ]
      # separate into metabolites, metabolites with adducts and isotopes
      grep_noiso_noadduct <- which(!grep("_", rownames(select_hmdb_df)))
      if (length(grep_noiso_noadduct) > 0) {
        select_metabolites <- select_hmdb_df[grep_noiso_noadduct, ]
      }
      grep_isotopes <- grep("iso", rownames(select_hmdb_df))
      if (length(grep_isotopes) > 0) {
        select_isotopes <- select_hmdb_df[grep_isotopes, ]
      }
      grep_adducts <- grep("_", rownames(select_hmdb_df[-grep_isotopes, ]))
      if (length(grep_adducts) > 0) {
        select_adducts <- select_hmdb_df[grep_adducts, ]
      }
      # take metabolite info first, then adducts. Isotope info in separate column.
      if (length(select_metabolites) > 0) {
        all_hmdb_names <- paste(select_metabolites[, "HMDB_name_all"], collapse = ";")
        all_hmdb_ids   <- paste(select_metabolites[, "HMDB_ID_all"], collapse = ";")
        sec_hmdb_ids   <- paste(select_metabolites[, "sec_HMDB_ID"], collapse = ";")
        theor_mz       <- select_metabolites[, column_label][1]
      }
      if (length(select_adducts) > 0) {
        all_hmdb_names <- paste(all_hmdb_names, paste(select_adducts[, "HMDB_name_all"], collapse =  ";"), collapse = ";")
        all_hmdb_ids   <- paste(all_hmdb_ids, paste(select_adducts[, "HMDB_ID_all"], collapse = ";"), collapse = ";")
        theor_mz <- select_adducts[, column_label][1]
      }
      if (length(select_isotopes) > 0) {
        all_isotope_names <- paste(select_isotopes[, "CompoundName"], collapse = ";")
      }
      assigned_hmdb[row_number, "assi_HMDB"] <- strsplit(all_hmdb_names, ";")[[1]][1]
      assigned_hmdb[row_number, "all_hmdb_names"] <- all_hmdb_names
      assigned_hmdb[row_number, "iso_HMDB"] <- all_isotope_names
      assigned_hmdb[row_number, "HMDB_code"] <- strsplit(all_hmdb_ids, ";")[[1]][1]
      assigned_hmdb[row_number, "all_hmdb_ids"] <- all_hmdb_ids
      assigned_hmdb[row_number, "sec_hmdb_ids"] <- sec_hmdb_ids
      assigned_hmdb[row_number, "theormz_HMDB"] <- theor_mz
    }
  }
  # combine all information
  peakgrouplist_identified <- cbind(ints_sorted[, 1:3], nrsamples = nrsamples, 
                                    ints_sorted[, 4:ncol(ints_sorted)], assigned_hmdb)
 
  return(peakgrouplist_identified)
}

