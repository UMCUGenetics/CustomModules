# unused variables will be removed: metab_list_alarm
prepare_data <- function(metab_list_all, zscore_patients_local) {
  #' Combine patient Z-scores with metabolite info
  #'
  #' @param metab_list_all: list of dataframes with metabolite information for different stofgroepen (list)
  #' @param zscore_patients_local: dataframe with metabolite Z-scores for all patient
  #'
  #' @return: list of dataframes for each stofgroep with data for each metabolite and patient/control per row

  # remove "_Zscore" from column (patient) names
  colnames(zscore_patients_local) <- gsub("_Zscore", "", colnames(zscore_patients_local))
  # put data into pages, max 20 violin plots per page in PDF
  metab_interest_sorted <- list()
  metab_category <- c()
  for (metab_class_index in 1:length(metab_list_all)) {
    metab_class <- names(metab_list_all)[metab_class_index]
    metab_list <- metab_list_all[[metab_class_index]]
    if (ncol(metab_list) > 2) {
      # third column are the alarm values, so reduce the data frame to 2 columns and save list
      metab_list_alarm <- metab_list
      metab_list <- metab_list[, c(1, 2)]
    }
    # make sure that all HMDB_names have 45 characters
    for (metab_index in 1:length(metab_list$HMDB_name)) {
      if (is.character(metab_list$HMDB_name[metab_index])) {
        hmdb_name_separated <- strsplit(metab_list$HMDB_name[metab_index], "")[[1]]
      } else {
        hmdb_name_separated <- "strspliterror"
      }
      if (length(hmdb_name_separated) <= 45) {
        hmdb_name_separated <- c(hmdb_name_separated, rep(" ", 45 - length(hmdb_name_separated)))
      } else {
        hmdb_name_separated <- c(hmdb_name_separated[1:42], "...")
      }
      metab_list$HMDB_name[metab_index] <- paste0(hmdb_name_separated, collapse = "")
    }
    # find metabolites and ratios in data frame zscore_patients_local
    metab_interest <- inner_join(metab_list, zscore_patients_local[-2], by = "HMDB_code")
    # remove column "HMDB_code"
    metab_interest <- metab_interest[, -which(colnames(metab_interest) == "HMDB_code")]
    # put the data frame in long format
    metab_interest_melt <- reshape2::melt(metab_interest, id.vars = "HMDB_name")
    # sort on metabolite names (HMDB_name)
    sort_order <- order(metab_interest_melt$HMDB_name)
    metab_interest_sorted[[metab_class_index]] <- metab_interest_melt[sort_order, ]
    metab_category <- c(metab_category, metab_class)
  }
  names(metab_interest_sorted) <- metab_category

  return(metab_interest_sorted)
}
