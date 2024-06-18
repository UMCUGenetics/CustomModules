is_diagnostic_patient <- function(patient_column){
  # Check for Diagnostics patients with correct patientnumber (e.g. starting with "P2024M")
  diagnostic_patients <- grepl("^P[0-9]{4}M",patient_column)
  
  return(diagnostic_patients)
}