# Functions for GenerateExcel

#' Get the indices of the control columns and a dataframe of intensities of the controls
#'
#' @param peakgroup_list: Dataframe with intensities for all samples (matrix)
#' @param label: Label used by grep to get the correct column (string)
#'
#' @returns: List with 2 items:
#'   col_idx: vector with indices of the control columns
#'   df_intensities: dataframe with the intensities of the controls
get_intensities_cols <- function(peakgroup_list, label) {
  # Get the column indices
  col_idx <- grep(label, colnames(peakgroup_list), fixed = TRUE)
  # Remove Z-score columns
  col_idx <- col_idx[!grepl("_Zscore", colnames(peakgroup_list)[col_idx], fixed = TRUE)]
  # Get intensity columns
  df_intensities <- peakgroup_list %>% select(all_of(col_idx))

  return(list(col_idx = col_idx, df_intensities = df_intensities))
}

#' Calculate the Z-scores with different methods for excluding outliers in controls
#'
#' @param peakgroup_list: Dataframe with intensities for all samples (matrix)
#' @param zscore_type: Method for excluding controls (string)
#' @param stat_filter: Either percentage or outlier threshold used for excluding controls (integer)
#' @param control_label: Part of name of all control samples (string)
#' @param case_label: Part of name of all patient samples (string)
#'
#' @returns: peakgroup_list_zscores: same dataframe as the input with added Z-score columns (matrix)
calculate_zscores <- function(peakgroup_list, zscore_type, stat_filter, control_label = "C", case_label = "P") {
  # Initialize
  peakgroup_list$avg_ctrls <- 0
  peakgroup_list$sd_ctrls <- 0
  
  # Get columns indices with intensities
  control_col_idx <- grep(control_label, colnames(peakgroup_list), fixed = TRUE)
  patient_col_idx <-  grep(case_label, colnames(peakgroup_list), fixed = TRUE)
  intensity_col_ids <- c(control_col_idx, patient_col_idx)
  peakgroup_list$nr_ctrls <- length(control_col_idx)
  
  # calculate mean and standard deviation of controls
  if (zscore_type == "_Zscore") {
    # Using all controls
    peakgroup_list$avg_ctrls <- apply(peakgroup_list[, control_col_idx], 1, function(x) mean(as.numeric(x), na.rm = TRUE))
    peakgroup_list$sd_ctrls <- apply(peakgroup_list[, control_col_idx], 1, function(x) sd(as.numeric(x), na.rm = TRUE))
  } else {
    if (length(control_col_idx) >= 3) {
      for (metabolite_index in seq_len(nrow(peakgroup_list))) {
        if (zscore_type == "_RobustZscore") {
          # remove outlier controls by using robust scaler
          peakgroup_list$avg_ctrls[metabolite_index] <- mean(robust_scaler(
            peakgroup_list[metabolite_index, control_col_idx],
            control_col_idx, stat_filter
          ))
          peakgroup_list$sd_ctrls[metabolite_index] <- sd(robust_scaler(
            peakgroup_list[metabolite_index, control_col_idx],
            control_col_idx, stat_filter
          ))
        } else {
          # remove outlier controls by using grubbs test
          intensities_without_outliers <- remove_outliers_grubbs(
            as.numeric(peakgroup_list[metabolite_index, control_col_idx]),
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
  peakgroup_list_zscores <- apply(peakgroup_list[, intensity_col_ids, drop = FALSE], 2, function(col) {
    (as.numeric(col) - peakgroup_list$avg_ctrls) / peakgroup_list$sd_ctrls
  })
  colnames(peakgroup_list_zscores) <- paste0(colnames(peakgroup_list)[intensity_col_ids], zscore_type)
  colnames(peakgroup_list_zscores) <- gsub("_OutlierRemovedZscore", "_Zscore", colnames(peakgroup_list_zscores))
  peakgroup_list_zscores <- cbind(peakgroup_list, peakgroup_list_zscores)
  
  return(peakgroup_list_zscores)
}

#' Robust scaler: remove outlier values in controls
#'
#' @param control_intensities: Intensities for control samples (vector of float)
#' @param control_col_ids: Column names for control samples (vector of string)
#' @param perc: Percentage of outliers which will be removed from controls (float)
#'
#' @returns trimmed_control_intensities: Intensities trimmed for outliers (vector of float)
robust_scaler <- function(control_intensities, control_col_ids, perc = 5) {
  # Determine how many values will be removed, based on percentage
  nr_to_remove <- ceiling(length(control_col_ids) * perc / 100)
  # Sort intensities
  sorted_control_intensities <- sort(as.numeric(control_intensities))
  # Remove highest and lowest intensities
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
  # Calculate Z-scores for all intensities
  mean_permetabolite <- mean(as.numeric(control_intensities))
  stdev_permetabolite <- sd(as.numeric(control_intensities))
  zscores_permetabolite <- (control_intensities - mean_permetabolite) / stdev_permetabolite
  # Remove intensities with a Z-score greater than outlier_threshold
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

#' Get the indices for intensity columns in a matrix
#'
#' @param peakgroup_list: Dataframe with intensities and possibly Z-scores for all samples (matrix)
#' @param control_label: part of name of all control samples (string)
#' @param case_label: part of name of all patient samples (string)
#'
#' @returns intensity_indices: indices of columns with intensities (vector of integers)
get_intensity_col_index <- function(peakgroup_list, control_label = "C", case_label = "P") {
  # Remove Zscore columns first
  if (any(grep("_Zscore", colnames(peakgroup_list)))) {
    peakgroup_list <- peakgroup_list[, -grep("_Zscore", colnames(peakgroup_list))]
  }
  # Get indices for controls
  control_col_ids <- grep(control_label, colnames(peakgroup_list), fixed = TRUE)
  # Get indices for patients
  patient_col_ids <- grep(case_label, colnames(peakgroup_list), fixed = TRUE)
  # Combine
  intensity_indices <- c(control_col_ids, patient_col_ids)

  return(intensity_indices)
}

#' Create an Excel workbook
#'
#' @param peakgroup_list: Dataframe with intensities for all samples (matrix)
#' @param intensity_col_ids: Indices for intensity columns (vector of integers)
#' @param prefix: Prefix to insert before file name (string)
#' @param z_score: Boolean number indicating whether Z-scores should be calculated (integer)
#' @param project: Name of dataset (string)
create_excel_output <- function(peakgroup_list, intensity_col_ids, prefix = "", z_score = 0, project) {
  # Set up
  if (prefix == "Drugs_") {
    plotdir <- "plots_drugs"
    sheetname <- "AllDrugs"
    wb_intensities <- openxlsx::createWorkbook("PeakGroups")
  } else if (prefix == "Helix_") {
    plotdir <- "plots_helix"
    sheetname <- "HelixSelection"
    wb_intensities <- openxlsx::createWorkbook("HelixOverview")
  } else {
    plotdir <- "plots"
    sheetname <- "BiologicallyRelevant"
    wb_intensities <- openxlsx::createWorkbook("AdductSums")
  }
  
  # Make sure the intensity columns are numeric
  peakgroup_list[intensity_col_ids] <- lapply(peakgroup_list[intensity_col_ids], as.numeric)
  # Round Z-scores to three decimals
  zscore_cols <- grep("_Zscore", colnames(peakgroup_list))
  peakgroup_list[zscore_cols] <- round(peakgroup_list[zscore_cols], 3)

  openxlsx::addWorksheet(wb_intensities, sheetname)

  # Add Z-scores and create plots
  if (z_score == 1) {
    dir.create(plotdir, showWarnings = FALSE)
    row_helix <- 2 # start on row 2 because of header
    # Get intensity columns
    intensities_df <- as.data.frame(peakgroup_list[ , intensity_col_ids])
    sample_names <- colnames(intensities_df)
    intensities_df$HMDB_key <- rownames(peakgroup_list)
    # Add a column for plots
    peakgroup_list <- cbind(plots = NA, peakgroup_list)

    for (row_index in seq_len(nrow(intensities_df))) {
      # Get HMDB ID
      hmdb_id <- intensities_df %>%
        slice(row_index) %>%
        pull(HMDB_key)

      # Transform dataframe to long format
      intensities_df_long <- intensities_df_to_long_format(intensities_df, row_index)

      # Set plot width to 40 times the number of samples
      plot_width <- length(unique(intensities_df_long$Samples)) * 40
      col_width <- plot_width * 2

      # Remove postive controls and SST mix samples, (e.g. P1001, P1002, P1003, P1005)
      intensities_df_long <- intensities_df_long %>% filter(!grepl("^P[0-9]{4}$", Samples))

      start_row_index <- row_index + 1
      save_plot_to_excel_workbook(
        wb_intensities,
        sheetname,
        intensities_df_long,
        paste0(plotdir, "/plot_"),
        hmdb_id,
        plot_width,
        col_width,
        start_row_index
      )
    }
    wb_intensities <- set_row_height_col_width_wb(
      wb_intensities,
      sheetname,
      nrow(peakgroup_list),
      ncol(peakgroup_list),
      col_width,
      plots_present = TRUE
    )

    # Reorder outlist for Excel file
    peakgroup_list <- peakgroup_list %>%
      relocate(c(HMDB_code, HMDB_name_all, avg_ctrls, sd_ctrls), .after = plots) %>%
      relocate(all_of(grep("_Zscore", colnames(peakgroup_list))), .after = sd_ctrls) %>%
      relocate(all_of(sample_names), .after = last_col())
    if (prefix == "Helix_") {
      colnames(peakgroup_list) <- gsub("H_Name", "Name", colnames(peakgroup_list))
      # Remove intensity columns
      peakgroup_list <- peakgroup_list %>% select(-all_of(sample_names))
    }
  } else if (z_score == 0) {
    save(peakgroup_list, file = "outlist.RData")
    if (!any(grepl("HMDB_code", colnames(peakgroup_list)))) {
      peakgroup_list$HMDB_code <- rownames(peakgroup_list)
    }
    wb_intensities <- set_row_height_col_width_wb(
      wb_intensities,
      sheetname,
      nrow(peakgroup_list),
      ncol(peakgroup_list),
      plot_width = NULL,
      plots_present = FALSE
    )
    peakgroup_list <- peakgroup_list %>%
      relocate(c(HMDB_name, HMDB_name_all, HMDB_code, HMDB_ID_all, sec_HMDB_ID))
    colnames(peakgroup_list) <- gsub("H_Name", "Name", colnames(peakgroup_list))
  }

  # Write Excel file
  openxlsx::writeData(wb_intensities, sheet = 1, peakgroup_list, startCol = 1)
  openxlsx::saveWorkbook(wb_intensities, paste0(prefix, project, ".xlsx"), overwrite = TRUE)
  rm(wb_intensities)
  unlink(plotdir, recursive = TRUE)
}

