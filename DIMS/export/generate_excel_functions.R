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
  # remove Z-score columns
  col_idx <- col_idx[!grepl("_Zscore", colnames(outlist)[col_idx], fixed = TRUE)]
  df_intensities <- as.data.frame(outlist[, col_idx])
  colnames(df_intensities) <- colnames(outlist)[col_idx]
  return(list(col_idx = col_idx, df_intensities = df_intensities))
}

#' Calculate the Z-scores with different methods for excluding outliers in controls
#'
#' @param peakgroup_list: Dataframe with intensities for all samples (matrix)
#' @param zscore_type: Method for excluding controls (string)
#' @param stat_filter: Either percentage or outlier threshold used for excluding controls (integer)
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
    # using all controls
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
  start_index <- nr_to_remove + 1
  end_index <- length(sorted_control_intensities) - nr_to_remove
  trimmed_control_intensities <- sorted_control_intensities[start_index:end_index]
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

#' Transform a dataframe with intensities to long format
#'
#' Get intensities of controls and patient for the selected metabolite,
#' pivot to long format, arrange Samples nummerically, change Sample names, get group size and
#' set Intensities to numeric.
#'
#' @param intensities_df: a dataframe with HMDB_key column and intensities for all samples
#'
#' @returns intensities_df_long: a dataframe with on each row a sample and their intensity
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
#' @param intensities_df_long: a dataframe with on each row a sample and their intensity
#' @param hmdb_id: HMDB ID of the selected metabolite
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

#' Make and save a boxplot of intensities to an Excel workbook
#'
#' For the Helix Excel the positive controls and SST mix samples are removed.
#'
#' @param excel_workbook: an openxlsx Workbook object
#' @param sheetname: a string containing the sheetname where the plots are to be placed
#' @param intensities_df: a dataframe containing intensities for controls and patients of a specific HMDB ID
#' @param file_path: a string containing the filepath for the png
#' @param hmdb_id: a string containing the HMDB ID that the intensities_df contains data for
#' @param plot_width: an integer containing the plot width for the png
#' @param col_width: an integer containing the width of the column that has the plots
#' @param start_row_index: an integer containing the index of the row where the plot has to be placed
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
  # remove Zscore columns first
  if (any(grep("_Zscore", colnames(peakgroup_list)))) {
    peakgroup_list <- peakgroup_list[, -grep("_Zscore", colnames(peakgroup_list))]
  }
  # get indices for controls
  control_col_ids <- grep(control_label, colnames(peakgroup_list), fixed = TRUE)
  # get indices for patients
  patient_col_ids <- grep(case_label, colnames(peakgroup_list), fixed = TRUE)
  # combine
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
  # set up
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

  openxlsx::addWorksheet(wb_intensities, sheetname)

  # Add Z-scores and create plots
  if (z_score == 1) {
    dir.create(plotdir, showWarnings = FALSE)
    row_helix <- 2 # start on row 2 because of header
    # get intensity columns
    intensities_df <- as.data.frame(peakgroup_list[ , intensity_col_ids])
    sample_names <- colnames(intensities_df)
    intensities_df$HMDB_key <- rownames(peakgroup_list)
    # add a column for plots
    peakgroup_list <- cbind(plots = NA, peakgroup_list)

    for (row_index in seq_len(nrow(intensities_df))) {
      # get HMDB ID
      hmdb_id <- intensities_df %>%
        slice(row_index) %>%
        pull(HMDB_key)

      # Transform dataframe to long format
      intensities_df_long <- intensities_df_to_long_format(intensities_df, row_index)

      # set plot width to 40 times the number of samples
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

    # reorder outlist for Excel file
    peakgroup_list <- peakgroup_list %>%
      relocate(c(HMDB_code, HMDB_name_all, avg_ctrls, sd_ctrls), .after = plots) %>%
      relocate(all_of(grep("_Zscore", colnames(peakgroup_list))), .after = sd_ctrls) %>%
      relocate(all_of(sample_names), .after = last_col())
    if (prefix == "Helix_") {
      colnames(peakgroup_list) <- gsub("H_Name", "Name", colnames(peakgroup_list))
      # remove intensity columns
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

  # write Excel file
  openxlsx::writeData(wb_intensities, sheet = 1, peakgroup_list, startCol = 1)
  openxlsx::saveWorkbook(wb_intensities, paste0(prefix, project, ".xlsx"), overwrite = TRUE)
  rm(wb_intensities)
  unlink(plotdir, recursive = TRUE)
}

