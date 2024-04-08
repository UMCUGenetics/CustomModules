prepare_alarmvalues <- function(pt_name, dims_helix_table) {
  # extract data for patient of interest (pt_name)
  pt_metabs_helix <- dims_helix_table %>% filter(Patient == pt_name)
  pt_metabs_helix$Z_score <- round(pt_metabs_helix$Z_score, 2)
  
  # Make empty dataframes for metabolites above or below alarmvalues
  pt_list_high <- data.frame(HMDB_name = character(), Z_score = numeric())
  pt_list_low <- data.frame(HMDB_name = character(), Z_score = numeric())
  
  # Loop over individual metabolites
  for (metab in unique(pt_metabs_helix$HMDB_name)){
    # Get data for individual metabolite
    pt_metab <- pt_metabs_helix %>% filter(HMDB_name == metab)
    # print(pt_metab)
    
    # Check if zscore is positive of negative
    if(pt_metab$Z_score > 0) {
      # Get specific alarmvalue for metabolite
      high_zscore_cutoff_metab <- pt_metabs_helix %>% filter(HMDB_name == metab) %>% pull(high_zscore)
      
      # If zscore is above the alarmvalue, add to pt_list_high table
      if(pt_metab$Z_score > high_zscore_cutoff_metab) {
        pt_metab_high <- pt_metab %>% select(HMDB_name, Z_score)
        pt_list_high <- rbind(pt_list_high, pt_metab_high)
      }
    } else {
      # Get specific alarmvalue for metabolite
      low_zscore_cutoff_metab <- pt_metabs_helix %>% filter(HMDB_name == metab) %>% pull(low_zscore)
      
      # If zscore is below the alarmvalue, add to pt_list_low table
      if(pt_metab$Z_score < low_zscore_cutoff_metab) {
        pt_metab_low <- pt_metab %>% select(HMDB_name, Z_score)
        pt_list_low <- rbind(pt_list_low, pt_metab_low)
      }
    }
  }
  
  # sort tables on zscore
  pt_list_high <- pt_list_high %>% arrange(desc(Z_score))
  pt_list_low <- pt_list_low %>% arrange(Z_score)
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
