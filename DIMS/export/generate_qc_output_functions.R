# functions for checking the quality of the dataset (internal standards, positive controls and SST sample)

#' Check the number of control samples for each metabolite and report findings in a txt file
#'
#' @param peakgroup_list: Dataframe with intensities and Z-scores for all samples (matrix)
#' @param min_num_controls: Minimum number of controls (integer)
#' @param file_name: Output file name (string)
check_number_of_controls <- function(peakgroup_list, min_num_controls, file_name) {
  peakgroups_below_nr_control <- peakgroup_list %>%
    filter(nr_ctrls < min_num_controls) %>%
    select(HMDB_name, HMDB_code)
  if (nrow(peakgroups_below_nr_control) == 0) {
    writeLines(paste0("All metabolites have", min_num_controls, "or more controls."), file_name)
  } else {
    writeLines(paste0("These metabolites have less than", min_num_controls, "controls. \n"), file_name)
    write.table(peakgroups_below_nr_control, file_name, append = TRUE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
}

#' Get the internal standards data
#'
#' @param internal_stand_df: dataframe with all data for all the internal standards (matrix)
#' @param scanmode: positive, negative or summed (string)
#' @param is_subset_filter: filter/threshold for outlier control removal (float)
#' @param dims_matrix: matrix used, e.g. Plasma, Research, etc. (string)
#' @param rundate: date of pipeline run (Date object)
#' @param project: project name (string)
#'
#' @returns internal_stand: dataframe with the intensities of the internal standards for all samples (matrix)
get_internal_standards <- function(internal_stand_df, scanmode, is_subset_filter, dims_matrix, rundate, project) {
  if (scanmode == "summed") {
    internal_stand <- internal_stand_df[c(names(is_subset_filter), "HMDB_code")]
    internal_stand$HMDB_name <- internal_stand_df$name
  } else {
    internal_stand <- as.data.frame(subset(is_subset_filter, rownames(is_subset_filter) %in% rownames(internal_stand_df)))
    internal_stand$HMDB_name <- internal_stand_df[
      match(row.names(internal_stand), internal_stand_df$HMDB_code, nomatch = NA),
      "name"
    ]
    internal_stand$HMDB_code <- row.names(internal_stand)
    internal_stand <- internal_stand %>% select(-c(HMDB_ID_all, sec_HMDB_ID, HMDB_name_all))
  }
  internal_stand <- reshape2::melt(internal_stand, id.vars = c("HMDB_code", "HMDB_name"))
  colnames(internal_stand) <- c("HMDB_code", "HMDB_name", "Sample", "Intensity")
  internal_stand$Matrix <- dims_matrix
  internal_stand$Rundate <- rundate
  internal_stand$Project <- project
  internal_stand$Intensity <- as.numeric(as.character(internal_stand$Intensity))

  return(internal_stand)
}

#' Generate and save internal standard plot
#'
#' @param plot_data: Dataframe with the data to be plotted (matrix)
#' @param plot_type: Type of plot (string)
#' @param plot_title: Title for the plot (string)
#' @param file_name: Name of the file (string)
#' @param plot_width: Width of the plot (int)
#' @param plot_height: Height of the plot (int)
#' @param hline_data: Dataframe with values for the minimal intensity line (matrix)
save_internal_standard_plot <- function(
    plot_data,
    plot_type,
    plot_title,
    file_name,
    plot_width,
    plot_height,
    hline_data = NULL) {

  # Check if plot_data contains data
  if (nrow(plot_data) == 0) {
    return()
  }

  if (plot_type == "barplot") {
    if (grepl("select", file_name)) {
      num_cols <- 2
    } else {
      num_cols <- NULL
    }

    plot <- ggplot2::ggplot(plot_data, aes(Sample_level, Intensity)) +
      ggplot2::geom_bar(aes(fill = HMDB_name), stat = "identity") +
      ggplot2::facet_wrap(~HMDB_name, scales = "free_y", ncol = num_cols) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::labs(x = "", y = "Intensity") +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6)
      )

    # Add a horizontal line if there is a minimum intensity given
    if (!is.null(hline_data)) {
      plot <- plot + ggplot2::geom_hline(aes(yintercept = int_line), subset(hline_data, HMDB_name %in% plot_data$HMDB_name))
    }
  } else if (plot_type == "lineplot") {
    plot <- ggplot2::ggplot(plot_data, aes(Sample_level, Intensity)) +
      ggplot2::geom_point(aes(col = HMDB_name)) +
      ggplot2::geom_line(aes(col = HMDB_name, group = HMDB_name)) +
      ggplot2::ggtitle(plot_title) +
      ggplot2::labs(x = "", y = "Intensity") +
      ggplot2::guides(
        shape = guide_legend(override.aes = list(size = 0.5)),
        color = guide_legend(override.aes = list(size = 0.5))
      ) +
      ggplot2::theme(
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.7, "line"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)
      )
  }

  ggplot2::ggsave(paste0("./plots/", file_name, ".png"),
    plot = plot, height = plot_height, width = plot_width, units = "in"
  )
}

#' Get the positive control data
#'
#' @param peakgroup_list: Dataframe with intensities for all samples (matrix)
#' @param sample_name: positive control sample name(s) (string)
#' @param hmdb_codes: HMDB codes of the positive control metabolites (vector)
#' @param hmdb_names: HMDB names of the positive control metabolites (vector)
#'
#' @returns: pos_ctrl_data: Dataframe with intensities and Z-scores of the positive control metabolites (matrix)
get_pos_ctrl_data <- function(peakgroup_list, sample_name, hmdb_codes, hmdb_names) {
  pos_ctrl_data <- peakgroup_list[hmdb_codes, c("HMDB_code", "name", sample_name)]
  pos_ctrl_data <- reshape2::melt(pos_ctrl_data, id.vars = c("HMDB_code", "name"))
  colnames(pos_ctrl_data) <- c("HMDB_code", "HMDB_name", "Sample", "Zscore")
  pos_ctrl_data$HMDB_name <- hmdb_names
  if ("Propionylglycine" %in% hmdb_names) {
    pos_ctrl_data$HMDB_code <- c("HMDB0000824", "HMDB0000783", "HMDB0000123")
  }
  return(pos_ctrl_data)
}

#' Round numbers to a set number of digits for numeric values
#'
#' @param df: Dataframe containing numeric values (matrix)
#' @param digits: Number of digits to round off to (integer)
#'
#' @return df: Dataframe with rounded numbers (matrix)
round_df <- function(df, digits) {
  numeric_columns <- sapply(df, mode) == "numeric"
  df[numeric_columns] <- round(df[numeric_columns], digits)
 
  return(df)
}

#' Get internal standard intensities
#'
#' @param is_data: Dataframe with intensities of internal standards (matrix)
#' @param int_cols: Indices of internal standard columns (vector of integers)
#' @param is_codes: Internal standard codes (vector of strings)
#'
#' @returns: is_intensities: Dataframe with intensities of internal standards (matrix)
get_is_intensities <- function(is_data, int_cols = NULL, is_codes = NULL) {
  if (is.null(is_codes)) {
    is_intensities <- is_data[, int_cols]
  } else {
    is_data <- as.data.frame(subset(is_data, rownames(is_data) %in% is_codes))
    is_intensities <- is_data %>% select(-c(HMDB_name, HMDB_ID_all, sec_HMDB_ID, HMDB_name_all))
  }
  is_intensities <- calc_coefficient_of_variation(is_intensities)
  is_intensities <- cbind(IS_name = is_data$HMDB_name, is_intensities)

  return(is_intensities)
}

#' Calculate coefficent of variation (cv) based on standard deviation (sd) and mean
#'
#' @param intensity_list: Intensities (matrix)
#'
#' @return intensity_list_with_cv: Intensities and cv, mean, sd (matrix)
calc_coefficient_of_variation <- function(intensity_list) {
  intensity_list <- as.data.frame(
    apply(intensity_list, 2, function(x) round(as.numeric(x), 0)),
    row.names = rownames(intensity_list)
  )
  sd_allsamples <- round(apply(intensity_list, 1, sd), 0)
  mean_allsamples <- round(apply(intensity_list, 1, mean), 0)
  cv_allsamples <- round(100 * sd_allsamples / mean_allsamples, 1)
  intensity_list_with_cv <- cbind(
    CV_perc = cv_allsamples,
    mean = mean_allsamples,
    sd = sd_allsamples,
    intensity_list
  )

  return(intensity_list_with_cv)
}

#' Check if all m/z values are present
#'
#' @param mzmed_pgrp_ident: All m/z values for a specific scanmode (vector of float)
#' @param scanmode: Scan mode, positive or negative (string)
#'
#' @return results_mz_missing: Either the missing mz values or a message that no values are missing (string)
check_missing_mz <- function(mzmed_pgrp_ident, scanmode) {
  # retrieve all unique m/z values in whole numbers and check if all are available
  mzmed_pgrp_ident <- unique(round(mzmed_pgrp_ident, digits = 0))
  # m/z range for a standard run = 70-600
  mz_range <- seq(70, 599, by = 1)
  mz_missing <- setdiff(mz_range, mzmed_pgrp_ident)
  # check if m/z are missing and make a txt file with information
  mz_missing_group <- cumsum(c(1, diff(mz_missing) != 1))
  if (length(mz_missing_group) > 1) {
    results_mz_missing <- c(paste0("Missing m/z values ", scanmode, " mode"))
    results_mz_missing <- c(results_mz_missing, by(mz_missing, mz_missing_group, identity))
  } else {
    results_mz_missing <- paste0(scanmode, " mode did not have missing mz values")
  }
  return(results_mz_missing)
}

#' Create a list of all internal standards with intensity below a threshold value
#'
#' @param is_selection_subset: Intensities for each internal standard in each sample (matrix)
#' @param thresholds: Threshold values for a given scan mode and matrix (vector of integers)
#' @param is_names: Names of internal standards for a given scan mode (vector of strings)
#' @param scanmode: Scan mode to include in output (string)
#'
#' @return is_below_threshold: All samples for which internal standard intensity is below threshold (matrix)
find_is_below_threshold <- function(is_selection_subset, thresholds, is_names, scanmode) {
  # for every line, check if intensity is below the appropriate threshold
  below_threshold_index <- c()
  for (line_index in seq_len(nrow(is_selection_subset))) {
    is_selected <- is_selection_subset$HMDB_name[line_index]
    thresh_selected <- thresholds[which(is_names == is_selected)]
    if (is_selection_subset$Intensity[line_index] < thresh_selected) {
      below_threshold_index <- c(below_threshold_index, line_index)
    }
  }
  is_below_threshold <- is_selection_subset[below_threshold_index, ]
  # add information on scan mode
  is_below_threshold <- cbind(is_below_threshold, scanmode = scanmode)
 
  return(is_below_threshold)
}

