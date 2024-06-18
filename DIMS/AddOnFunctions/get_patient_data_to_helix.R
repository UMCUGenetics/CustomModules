get_patient_data_to_helix <- function(metab_interest_sorted, metab_list_all){
  # Combine Z-scores of metab groups together
  df_all_metabs_zscores <- bind_rows(metab_interest_sorted)
  # Change columnnames 
  colnames(df_all_metabs_zscores) <- c("HMDB_name", "Patient", "Z_score")
  # Change Patient column to character instead of factor
  df_all_metabs_zscores$Patient <- as.character(df_all_metabs_zscores$Patient)
  
  # Delete whitespaces HMDB_name 
  df_all_metabs_zscores$HMDB_name <- str_trim(df_all_metabs_zscores$HMDB_name, "right")
  
  # Split HMDB_name column on "nitine;" for match dims_helix_table
  df_all_metabs_zscores$HMDB_name_split <- str_split_fixed(df_all_metabs_zscores$HMDB_name, "nitine;", 2)[,1]
  
  # Combine stofgroepen
  dims_helix_table <- bind_rows(metab_list_all)
  # Filter table for metabolites for Helix
  dims_helix_table <- dims_helix_table %>% filter(Helix == "ja")
  # Split HMDB_name column on "nitine;" for match df_all_metabs_zscores
  dims_helix_table$HMDB_name_split <- str_split_fixed(dims_helix_table$HMDB_name, "nitine;", 2)[,1]
  dims_helix_table <- dims_helix_table %>% select(HMDB_name_split, Helix_naam, high_zscore, low_zscore)
  
  # Filter DIMS results for metabolites for Helix
  df_metabs_helix <- df_all_metabs_zscores %>% filter(HMDB_name_split %in% dims_helix_table$HMDB_name_split)
  # Combine dims_helix_table and df_metabs_helix, adding Helix codes etc.
  df_metabs_helix <- df_metabs_helix %>% left_join(dims_helix_table, by = join_by(HMDB_name_split))
  
  df_metabs_helix <- df_metabs_helix %>% select(HMDB_name, Patient, Z_score, Helix_naam, high_zscore, low_zscore)
  
  return(df_metabs_helix)
}