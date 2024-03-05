prepare_toplist <- function(pt_name, zscore_patients_copy) {
  # set parameters for table
  top_highest <- 20
  top_lowest <- 10
  
  # extract data for patient of interest (pt_name)
  pt_list <- zscore_patients_copy[ , c(1,2, which(colnames(zscore_patients_copy) == pt_name))]
  # sort metabolites on Z-scores for this patient
  pt_list_sort <- sort(pt_list[ , 3], index.return=TRUE)
  pt_list_sorted <- pt_list[pt_list_sort$ix, ]
  # determine top highest and lowest Z-scores for this patient
  pt_list_sort <- sort(pt_list[ , 3], index.return=TRUE)
  pt_list_low  <- pt_list[pt_list_sort$ix[1:top_lowest], ]
  pt_list_high <- pt_list[pt_list_sort$ix[length(pt_list_sort$ix):(length(pt_list_sort$ix)-top_highest+1)], ]
  # round off Z-scores
  pt_list_low[ , 3]  <- round(as.numeric(pt_list_low[ , 3]), 2)
  pt_list_high[ , 3] <- round(as.numeric(pt_list_high[ , 3]), 2)
  # add lines for increased, decreased
  extra_line1 <- c("Increased", "", "") 
  extra_line2 <- c("Decreased", "", "")
  top_metab_pt <- rbind(extra_line1, pt_list_high, extra_line2, pt_list_low)
  # remove row names
  rownames(top_metab_pt) <- NULL
  
  # change column names for display
  colnames(top_metab_pt) <- c("HMDB_ID", "Metabolite", "Z-score")
  
  return(top_metab_pt)
}
