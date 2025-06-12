# Functions for GenerateExcel
get_intensities_cols <- function(outlist, label) {
  #' Get the indices of the control columns and a dataframe of intensities of the controls
  #'
  #' @param outlist: dataframe with intensities for all samples
  #' @param label: string used by grep to get the correct columns
  #'
  #' @returns: list with 2 items:
  #'                  col_idx: vector with indices of the control columns
  #'                  df_intensities: dataframe with the intensities of the controls
  col_idx <- grep(label, colnames(outlist), fixed = TRUE)
  df_intensities <- as.data.frame(outlist[, col_idx])
  colnames(df_intensities) <- colnames(outlist)[col_idx]
  return(list(col_idx = col_idx, df_intensities = df_intensities))
}

calculate_zscores <- function(outlist, zscore_type, control_cols, stat_filter, intensity_col_ids, startcol) {
  #' Calculate the Z-scores with different methods for excluding controls
  #'
  #' @param outlist: dataframe with intensities for all samples
  #' @param zscore_type: string with method for excluding controls
  #' @param control_cols: vector with indices of the control columns
  #' @param stat_filter: integer used for excluding controls, either percentage or outlier threshold
  #' @param intensity_col_ids: vector with indices of the samples for which to calculate Z-scores
  #' @param startcol: integer of the column from where to add the Z-score columns
  #'
  #' @returns: outlist: same dataframe as the input with added Z-score columns

  # Calculate mean and sd
  outlist$avg_ctrls <- 0
  outlist$sd_ctrls <- 0
  outlist$nr_ctrls <- length(control_cols)

  if (zscore_type == "_Zscore") {
    # Calculate mean and sd with all controls
    outlist$avg_ctrls <- apply(control_cols, 1, function(x) mean(as.numeric(x), na.rm = TRUE))
    outlist$sd_ctrls  <- apply(control_cols, 1, function(x) sd(as.numeric(x), na.rm = TRUE))
  } else {
    if (length(control_cols) > 3) {
      for (metabolite_index in seq_len(nrow(outlist))) {
        if (zscore_type == "_RobustZscore") {
          # Calculate mean and sd, remove outlier controls by using robust scaler
          outlist$avg_ctrls[metabolite_index] <- mean(robust_scaler(outlist[metabolite_index, control_cols],
                                                                    control_cols, stat_filter))
          outlist$sd_ctrls[metabolite_index]  <-   sd(robust_scaler(outlist[metabolite_index, control_cols],
                                                                    control_cols, stat_filter))
        } else {
          # Calculate mean, sd and number of remaining controls, remove outlier controls by using grubbs test
          intensities_without_outliers <- remove_outliers_grubbs(as.numeric(outlist[metabolite_index, control_cols]), stat_filter)
          outlist$avg_ctrls[metabolite_index] <- mean(intensities_without_outliers)
          outlist$sd_ctrls[metabolite_index]  <- sd(intensities_without_outliers)
          outlist$nr_ctrls[metabolite_index]  <- length(intensities_without_outliers)
        }
      }
    }
  }

  # Calculate Z-scores
  outlist_zscores <- apply(outlist[, intensity_col_ids, drop = FALSE], 2, function(col) {
    (as.numeric(col) - outlist$avg_ctrls) / outlist$sd_ctrls
  })
  outlist <- cbind(outlist, outlist_zscores)
  colnames(outlist)[startcol:ncol(outlist)] <- paste0(colnames(outlist)[intensity_col_ids], zscore_type)
  
  return(outlist)
}

robust_scaler <- function(control_intensities, control_col_ids, perc = 5) {
  #' Calculate robust scaler: Z-score based on controls without outliers
  #'
  #' @param control_intensities: Matrix with intensities for control samples
  #' @param control_col_ids: Vector with column names for control samples
  #' @param perc: Percentage of outliers which will be removed from controls (float)
  #'
  #' @return trimmed_control_intensities: Intensities trimmed for outliers
  nr_to_remove <- ceiling(length(control_col_ids) * perc / 100)
  sorted_control_intensities <- sort(as.numeric(control_intensities))
  trimmed_control_intensities <- sorted_control_intensities[(nr_to_remove + 1) :
                                                              (length(sorted_control_intensities) - nr_to_remove)]
  return(trimmed_control_intensities)
}

remove_outliers_grubbs <- function(control_intensities, outlier_threshold = 2) {
  #' Remove outliers per metabolite according to Grubb's test
  #'
  #' @param control_intensities: Vector with intensities for control samples
  #' @param outlier_threshold: Threshold for outliers which will be removed from controls (float)
  #'
  #' @return trimmed_control_intensities: Intensities trimmed for outliers
  mean_permetabolite <- mean(as.numeric(control_intensities))
  stdev_permetabolite <- sd(as.numeric(control_intensities))
  zscores_permetabolite <- (control_intensities - mean_permetabolite) / stdev_permetabolite
  # remove intensities with a zscore_permetabolite greater than outlier_threshold
  if (sum(zscores_permetabolite > outlier_threshold) > 0) {
    trimmed_control_intensities <- control_intensities[-which(zscores_permetabolite > outlier_threshold)]
  } else {
    trimmed_control_intensities <- control_intensities
  }
  return(trimmed_control_intensities)
}

save_to_rdata_and_txt <- function(df, file_name) {
  #' Save a dataframe to RData and txt
  #'
  #' @param df: dataframe
  #' @param file_name: string with the file name
  save(df, file = paste0(file_name, ".RData"))
  write.table(df, file = paste0(file_name, ".txt"), sep = "\t", row.names = FALSE)
}

set_row_height_col_width_wb <- function(wb, sheetname, num_rows_df, num_cols_df, plot_width, plots_present) {
  #' Change the row height and column width of the Excel
  #'
  #' @param wb: an openxlsx workbook (S4 object)
  #' @param sheetname: name of the workbook sheet (string)
  #' @param num_rows_df: number of rows in dataframe (int)
  #' @param num_col_df: number of columns in dataframe (int)
  #' @param plot_width: width of the plots to be added (int)
  #' @param plots_present: boolean if plots are added to the workbook (boolean)
  #'
  #' @returns wb: a workbook object with changed row heights and column widths
  if (plots_present) {
    openxlsx::setColWidths(wb, sheetname, cols = 1, widths = plot_width / 20)
    openxlsx::setRowHeights(wb, sheetname, rows = c(seq(2, num_rows_df + 1)), heights = 560 / 4)
    openxlsx::setColWidths(wb, sheetname, cols = c(seq(2, num_cols_df)), widths = 20)
  } else {
    openxlsx::setRowHeights(wb, sheetname, rows = c(seq_len(num_rows_df)), heights = 18)
    openxlsx::setColWidths(wb, sheetname, cols = c(seq_len(num_cols_df)), widths = 20)
  }
  return(wb)
}
