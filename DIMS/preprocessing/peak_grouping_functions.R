# functions for PeakGrouping
find_peak_groups <- function(outlist_sorted, mz_tolerance, sample_names) {
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
      select_peaks <- select_peaks %>% dplyr::group_by(samplenr) %>% dplyr::slice_min(order_by = abs_ppmdev, n = 1)
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

annotate_peak_groups <- function(ints_sorted, hmdb_add_iso) {
  #' annotate peak groups; assign metabolites (adducts, isotopes) with suitable mass from HMDB
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

