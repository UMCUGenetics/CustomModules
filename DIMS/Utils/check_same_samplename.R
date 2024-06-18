check_same_samplename <- function(int_col_name, zscore_col_name) {
  #' A check to see if intensity and Z-score columns match
  #'
  #' @param int_col_name: name of an intensity column (string)
  #' @param zscore_col_name: name of a Z-score column (string)
  #'
  #' @return: match or mismatch of the columns (boolean)

  paste0(int_col_name, "_Zscore") == zscore_col_name
}
