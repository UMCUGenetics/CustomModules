add_lab_id_and_onderzoeksnummer <- function(df_metabs_helix) {
  # Split patient number into labnummer and Onderzoeksnummer
  for (row in 1:nrow(df_metabs_helix)) {
    df_metabs_helix[row,"labnummer"] <- gsub("^P|\\.[0-9]*", "", df_metabs_helix[row,"Patient"])
    labnummer_split <- strsplit(as.character(df_metabs_helix[row, "labnummer"]), "M")[[1]]
    df_metabs_helix[row, "Onderzoeksnummer"] <- paste0("MB", labnummer_split[1], "/", labnummer_split[2])
  }
  
  return(df_metabs_helix)
}
