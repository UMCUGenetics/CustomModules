# Functions for GenerateExcel

#' Get the indices of the control columns and a dataframe of intensities of the controls
#'
#' @param peakgroup_list: Dataframe with intensities for all samples (matrix)
#' @param label: Label used by grep to get the correct column (string)
#'
#' @returns: List with 2 items:
#'                  col_idx: vector with indices of the control columns
#'                  df_intensities: dataframe with the intensities of the controls
get_intensities_cols <- function(peakgroup_list, label) {
  # get the column indices
  col_idx <- grep(label, colnames(peakgroup_list), fixed = TRUE)
  # get intensity columns
  df_intensities <- as.data.frame(peakgroup_list[, col_idx])
  colnames(df_intensities) <- colnames(peakgroup_list)[col_idx]

  return(list(col_idx = col_idx, df_intensities = df_intensities))
}

#' Calculate the Z-scores with different methods for excluding outliers in controls
#'
#' @param peakgroup_list: dataframe with intensities for all samples (matrix)
#' @param zscore_type: string with method for excluding controls
#' @param control_cols: vector with indices of the control columns
#' @param stat_filter: integer used for excluding controls, either percentage or outlier threshold
#' @param intensity_col_ids: vector with indices of the samples for which to calculate Z-scores
#'
#' @returns: peakgroup_list: same dataframe as the input with added Z-score columns (matrix)
calculate_zscores <- function(peakgroup_list, zscore_type, control_cols, stat_filter, intensity_col_ids) {

  # Calculate mean and sd
  peakgroup_list$avg_ctrls <- 0
  peakgroup_list$sd_ctrls <- 0
  peakgroup_list$nr_ctrls <- length(control_cols)

  if (zscore_type == "_Zscore") {
    # Calculate mean and sd with all controls
    peakgroup_list$avg_ctrls <- apply(control_cols, 1, function(x) mean(as.numeric(x), na.rm = TRUE))
    peakgroup_list$sd_ctrls <- apply(control_cols, 1, function(x) sd(as.numeric(x), na.rm = TRUE))
  } else {
    if (length(control_cols) > 3) {
      for (metabolite_index in seq_len(nrow(peakgroup_list))) {
        if (zscore_type == "_RobustZscore") {
          # Calculate mean and sd, remove outlier controls by using robust scaler
          peakgroup_list$avg_ctrls[metabolite_index] <- mean(robust_scaler(
            peakgroup_list[metabolite_index, control_cols],
            control_cols, stat_filter
          ))
          peakgroup_list$sd_ctrls[metabolite_index] <- sd(robust_scaler(
            peakgroup_list[metabolite_index, control_cols],
            control_cols, stat_filter
          ))
        } else {
          # Calculate mean, sd and number of remaining controls, remove outlier controls by using grubbs test
          intensities_without_outliers <- remove_outliers_grubbs(
            as.numeric(peakgroup_list[metabolite_index, control_cols]),
            stat_filter
          )
          peakgroup_list$avg_ctrls[metabolite_index] <- mean(intensities_without_outliers)
          peakgroup_list$sd_ctrls[metabolite_index] <- sd(intensities_without_outliers)
          peakgroup_list$nr_ctrls[metabolite_index] <- length(intensities_without_outliers)
        }
      }
    }
  }

  # Calculate Z-scores
  outlist_zscores <- apply(peakgroup_list[, intensity_col_ids, drop = FALSE], 2, function(col) {
    (as.numeric(col) - peakgroup_list$avg_ctrls) / peakgroup_list$sd_ctrls
  })
  colnames(peakgroup_list) <- paste0(colnames(peakgroup_list)[intensity_col_ids], zscore_type)
  peakgroup_list <- cbind(peakgroup_list, outlist_zscores)

  return(peakgroup_list)
}

#' Robust scaler: remove outlier values in controls
#'
#' @param control_intensities: Intensities for control samples (vector of float)
#' @param control_col_ids: Column names for control samples (vector of string)
#' @param perc: Percentage of outliers which will be removed from controls (float)
#'
#' @returns trimmed_control_intensities: Intensities trimmed for outliers (vector of float)
robust_scaler <- function(control_intensities, control_col_ids, perc = 5) {
  # determine how many values will be removed, based on percentage
  nr_to_remove <- ceiling(length(control_col_ids) * perc / 100)
  # sort intensities
  sorted_control_intensities <- sort(as.numeric(control_intensities))
  # remove highest and lowest intensities
  start_index <- nr_to_remove + 1
  end_index <- length(sorted_control_intensities) - nr_to_remove
  trimmed_control_intensities <- sorted_control_intensities[start_index:end_index]

  return(trimmed_control_intensities)
}

#' Remove outlier values using Grubb's test
#'
#' @param control_intensities: Intensities for control samples (vector of float)
#' @param outlier_threshold: Threshold for outliers to be removed from controls (float)
#'
#' @returns trimmed_control_intensities: Intensities trimmed for outliers (vector of float)
remove_outliers_grubbs <- function(control_intensities, outlier_threshold = 2) {
  # calculate Z-scores for all intensities
  mean_permetabolite <- mean(as.numeric(control_intensities))
  stdev_permetabolite <- sd(as.numeric(control_intensities))
  zscores_permetabolite <- (control_intensities - mean_permetabolite) / stdev_permetabolite
  # remove intensities with a Z-score greater than outlier_threshold
  if (sum(zscores_permetabolite > outlier_threshold) > 0) {
    trimmed_control_intensities <- control_intensities[-which(zscores_permetabolite > outlier_threshold)]
  } else {
    trimmed_control_intensities <- control_intensities
  }

  return(trimmed_control_intensities)
}

#' Save a dataframe to RData and txt
#'
#' @param df: Dataframe (matrix)
#' @param file_name: File name (string)
save_to_rdata_and_txt <- function(df, file_name) {
  save(df, file = paste0(file_name, ".RData"))
  write.table(df, file = paste0(file_name, ".txt"), sep = "\t", row.names = FALSE)
}

#' Set the row height and column width of the Excel
#'
#' @param wb: An openxlsx workbook (S4 object)
#' @param sheetname: Name of the workbook sheet (string)
#' @param num_rows_df: Number of rows in dataframe (integer)
#' @param num_col_df: Number of columns in dataframe (integer)
#' @param plot_width: Width of the plots to be added (integer)
#' @param plots_present: Parameter that indicates whether plots are added to the workbook (boolean)
#'
#' @returns wb: Workbook object with changed row heights and column widths
set_row_height_col_width_wb <- function(wb, sheetname, num_rows_df, num_cols_df, plot_width, plots_present) {
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

#' Transform a dataframe with intensities to long format
#'
#' Get intensities of controls and patient for the selected metabolite,
#' pivot to long format, arrange Samples nummerically, change Sample names, get group size and
#' set Intensities to numeric.
#'
#' @param intensities_df: Dataframe with HMDB_key column and intensities (matrix) 
#' @param row_index: Index of row (integer) 
#'
#' @returns intensities_df_long: Dataframe with on each row a sample and its intensity (matrix)
intensities_df_to_long_format <- function(intensities_df, row_index) {
  intensities_df_long <- intensities_df %>%
    slice(row_index) %>%
    select(-HMDB_key) %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "Samples", values_to = "Intensities") %>%
    arrange(nchar(Samples)) %>%
    mutate(
      Samples = gsub("\\..*", "", Samples),
      Samples = gsub("(C).*", "\\1", Samples),
      Intensities = as.numeric(Intensities),
      type = ifelse(Samples == "C", "Control", "Patients")
    ) %>%
    group_by(Samples) %>%
    mutate(group_size = n()) %>%
    ungroup()

  return(intensities_df_long)
}

#' Create a plot of intensities of samples for Excel
#' Use boxplot if group size is above 2, otherwise use a dash/line
#'
#' @param intensities_df_long: Dataframe with on each row a sample and its intensity (matrix)
#' @param hmdb_id: HMDB ID of the selected metabolite (string)
#'
#' @returns boxplot_object: ggplot2 object containing the plot of intensities
create_boxplot <- function(intensities_df_long, hmdb_id) {
  boxplot_object <- ggplot(intensities_df_long, aes(Samples, Intensities)) +
    geom_boxplot(data = subset(intensities_df_long, group_size > 2), aes(fill = type)) +
    geom_point(
      data = subset(intensities_df_long, group_size <= 2),
      shape = "-",
      size = 10,
      aes(colour = type, fill = type)
    ) +
    scale_fill_manual(values = c("Control" = "green", "Patients" = "#b20000")) +
    scale_color_manual(values = c("Control" = "black", "Patients" = "#b20000")) +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      panel.background = element_rect(fill = "white", colour = "black")
    ) +
    ggtitle(hmdb_id)

  return(boxplot_object)
}

#' Make and save a boxplot of intensities in png format and insert into an Excel workbook
#' For the Helix Excel, the positive controls and SST mix samples are removed.
#'
#' @param excel_workbook: An openxlsx Workbook object (workbook object)
#' @param sheetname: Name of the sheet where the plots are to be placed (string)
#' @param intensities_df: Dataframe containing intensities for controls and patients of a specific HMDB ID (matrix)
#' @param file_path: Filepath for the png (string)
#' @param hmdb_id: HMDB ID corresponding to intensities_df (string) 
#' @param plot_width: Plot width for the png (integer)
#' @param col_width: Width of the column that has the plots (integer)
#' @param start_row_index: Index of the row where the plot has to be placed (integer)
save_plot_to_excel_workbook <- function(excel_workbook,
                                        sheetname,
                                        intensities_df,
                                        file_path,
                                        hmdb_id,
                                        plot_width,
                                        col_width,
                                        start_row_index) {
  plot.new()
  tmp_png <- paste0(file_path, hmdb_id, ".png")
  png(filename = tmp_png, width = plot_width, height = 300)

  boxplot <- create_boxplot(intensities_df, hmdb_id)

  print(boxplot)
  dev.off()

  openxlsx::insertImage(
    excel_workbook,
    sheetname,
    tmp_png,
    startRow = start_row_index,
    startCol = 1,
    height = 560,
    width = col_width,
    units = "px"
  )

  return(excel_workbook)
}
