prepare_alarmvalues <- function(pt_name, metab_interest_sorted) {
  # set parameters for table
  high_zscore_cutoff <- 5
  low_zscore_cutoff <- -3
  
  # make table of all metabolites
  all_metab <- c()
  for (page_nr in 1:length(metab_interest_sorted)) {
    all_metab <- rbind(all_metab, metab_interest_sorted[[page_nr]])
  }
  # extract data for patient of interest (pt_name)
  pt_list <- all_metab[which(all_metab$variable==pt_name), ]
  # remove column with patient name
  pt_list <- pt_list[ , -2]
  # round off Z-scores
  pt_list$value <- round(as.numeric(pt_list$value), 2)
  
  # determine alarms for this patient: Z-score above 5 or below -3
  pt_list_high <- pt_list[pt_list$value > high_zscore_cutoff, ]
  pt_list_low <- pt_list[pt_list$value < low_zscore_cutoff, ]
  # add lines for increased, decreased
  extra_line1 <- c("Increased", "") 
  extra_line2 <- c("Decreased", "")
  # combine the two lists
  top_metab_pt <- rbind(extra_line1, pt_list_high, extra_line2, pt_list_low)
  # remove row names
  rownames(top_metab_pt) <- NULL
  # change column names for display
  colnames(top_metab_pt) <- c("Metabolite", "Z-score")

  return(top_metab_pt)
}
