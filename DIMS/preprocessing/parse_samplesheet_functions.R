# function for parse_samplesheet
generate_repl_pattern <- function(sample_sheet) {
  #' Generate replication pattern list based on information in sample_sheet
  #'
  #' @param sample_names: vector of sample names (vector of strings)
  #' @param sample_sheet: matrix of file names and sample names
  #'
  #' @return ints_sorted: list of sample names with corresponding file names (technical replicates)

  # get the right columns from the samplesheet
  file_name_col <- grep("File_Name|File Name", colnames(sample_sheet))
  sample_name_col <- grep("Sample_Name|Sample Name", colnames(sample_sheet))
  # get the unique sample names from the samplesheet
  sample_names <- sort(unique(trimws(as.vector(unlist(sample_sheet[sample_name_col])))))
  # remove all characters from sample_names which are not letters, numbers, hyphens and periods
  sample_names <- gsub("[^-.[:alnum:]]", "_", sample_names)

  # create replication pattern (which technical replicates belong to which sample)
  repl_pattern <- c()
  for (sample_group in sample_names) {
    file_indices <- which(sample_sheet[, sample_name_col] == sample_group)
    file_names <- sample_sheet[file_indices, file_name_col]
    repl_pattern <- c(repl_pattern, list(file_names))
  }
  names(repl_pattern) <- sample_names

  return(repl_pattern)
}

