# Function for parsing samplesheet

#' Generate replication pattern list based on information in sample_sheet
#'
#' @param sample_sheet: File names and sample names (matrix)
#'
#' @return ints_sorted: Sample names with corresponding file names (list of strings))
generate_repl_pattern <- function(sample_sheet) {
  # Get the file name and sample name columns from the samplesheet
  file_name_col <- grep("File_Name|File Name", colnames(sample_sheet), ignore.case = TRUE)
  sample_name_col <- grep("Sample_Name|Sample Name", colnames(sample_sheet), ignore.case = TRUE)
  # Get the unique sample names from the samplesheet
    sample_names <- sample_sheet[sample_name_col] |>
    unlist() |>
    as.vector() |>
    trimws() |>
    unique() |>
    sort()
  # Remove all characters from sample_names which are not letters, numbers, hyphens and periods
  sample_names <- gsub("[^-.[:alnum:]]", "_", sample_names)

  # Create replication pattern (which technical replicates belong to which sample)
  repl_pattern <- split(
    sample_sheet[[file_name_col]],
    sample_sheet[[sample_name_col]]
  )

  return(repl_pattern)
}
