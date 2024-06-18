getPatients <- function(peaklist){
  patients=colnames(peaklist)[grep("P", colnames(peaklist), fixed = TRUE)]
  patients=unique(as.vector(unlist(lapply(strsplit(patients, ".", fixed = TRUE), function(x) x[1]))))
  # ToDo: If 2 P's in sample names!!!!!!!!!!!!!
  # patients=sort(as.numeric(unique(as.vector(unlist(lapply(strsplit(patients, "_P", fixed = TRUE), function(x) x[2]))))))
  patients=sort(as.numeric(unique(as.vector(unlist(lapply(strsplit(patients, "P", fixed = TRUE), function(x) x[2]))))))
  
  return(patients)
}
