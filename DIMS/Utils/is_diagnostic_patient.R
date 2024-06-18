is_diagnostic_patient <- function(patient_column) {
  #' Check for Diagnostics patients with correct patient number (e.g. starting with "P2024M")
  #'
  #' @param patient_column: a column from dataframe with IDs (character vector)
  #'
  #' @return: a logical vector with TRUE or FALSE for each element (vector)
  
  diagnostic_patients <- grepl("^P[0-9]{4}M", patient_column)
  
  return(diagnostic_patients)
}
