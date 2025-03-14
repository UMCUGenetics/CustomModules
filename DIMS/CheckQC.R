library("ggplot2")
library("reshape2")
library("openxlsx")
suppressMessages(library("dplyr"))

# set the number of digits for floats
options(digits = 16)

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

init_file <- cmd_args[1]
project <- cmd_args[2]
dims_matrix <- cmd_args[3]
z_score <- cmd_args[4]
sst_components_file <- cmd_args[5]
export_scripts_dir <- cmd_args[6]

outdir <- "./"

# load in function scripts
source(paste0(export_scripts_dir, "check_qc_functions.R"))

# load init files
load(init_file)
# load outlist from GenerateExcel
load("outlist.RData")
# load combined adducts for each scanmodus
load("AdductSums_positive.RData")
load("AdductSums_negative.RData")

# get current date
rundate <- Sys.Date()

# create a directory for plots in project directory
dir.create(paste0(outdir, "/plots"), showWarnings = FALSE)

control_label <- "C"

#### INTERNAL STANDARDS ####
internal_stand_list <- outlist[grep("Internal standard", outlist[, "relevance"], fixed = TRUE), ]
internal_stand_codes <- rownames(internal_stand_list)

# check if there is data present for all the samples that the pipeline started with,
# if not write sample name to a log file.
sample_names_nodata <- setdiff(names(repl_pattern), names(internal_stand_list))
if (!is.null(sample_names_nodata)) {
  write.table(sample_names_nodata, file = paste(outdir, "sample_names_nodata.txt", sep = "/"),
              rowNames = FALSE, colNames = FALSE, quote = FALSE)
  for (sample_name in sample_names_nodata) {
    repl_pattern[[sample_name]] <- NULL
  }
}

# Only continue with patients (columns) that are in both pos and neg, so patients that are in both
samples_both_modes <- intersect(colnames(outlist_tot_neg), colnames(outlist_tot_pos))
outlist_tot_neg <- outlist_tot_neg[, samples_both_modes]
outlist_tot_pos <- outlist_tot_pos[, samples_both_modes]

# Retrieve IS summed adducts
internal_stand_summed <- get_internal_standards(internal_stand_list, "summed", repl_pattern, dims_matrix, rundate, project)
# Retrieve IS positive mode
internal_stand_pos <- get_internal_standards(internal_stand_list, "pos", outlist_tot_pos, dims_matrix, rundate, project)
# Retrieve IS negative mode
internal_stand_neg <- get_internal_standards(internal_stand_list, "neg", outlist_tot_neg, dims_matrix, rundate, project)

# Save results
save(internal_stand_pos, internal_stand_neg, internal_stand_summed, file = paste0(outdir, "/", project, "_IS_results.RData"))

# number of samples, for plotting length and width
sample_count <- length(repl_pattern)

# change the order of the x-axis summed plots to a natural sorted one
sample_naturalorder <- unique(as.character(internal_stand_summed$Sample))
sample_naturalorder <- stringr::str_sort(sample_naturalorder, numeric = TRUE)
internal_stand_summed$Sample_level <- factor(internal_stand_summed$Sample, levels = c(sample_naturalorder))
internal_stand_pos$Sample_level <- factor(internal_stand_pos$Sample, levels = c(sample_naturalorder))
internal_stand_neg$Sample_level <- factor(internal_stand_neg$Sample, levels = c(sample_naturalorder))

## bar plots with all IS
plot_width <- 9 + 0.35 * sample_count
plot_height <- plot_width / 2.5

save_internal_standard_plot(internal_stand_neg, "barplot", "Interne Standaard (Neg)", outdir, "IS_bar_all_neg", plot_width, plot_height)
save_internal_standard_plot(internal_stand_pos, "barplot", "Interne Standaard (Pos)", outdir, "IS_bar_all_pos", plot_width, plot_height)
save_internal_standard_plot(internal_stand_summed, "barplot", "Interne Standaard (Summed)", outdir, "IS_bar_all_sum", plot_width, plot_height)

## Line plots with all IS
plot_width <- 8 + 0.2 * sample_count
plot_height <- plot_width / 2.5

save_internal_standard_plot(internal_stand_neg, "lineplot", "Interne Standaard (Neg)", outdir, "IS_line_all_neg", plot_width, plot_height)
save_internal_standard_plot(internal_stand_pos, "lineplot", "Interne Standaard (Pos)", outdir, "IS_line_all_pos", plot_width, plot_height)
save_internal_standard_plot(internal_stand_summed, "lineplot", "Interne Standaard (Sum)", outdir, "IS_line_all_sum", plot_width, plot_height)

## bar plots with a selection of IS
internal_stand_neg_selection <- c("2H2-Ornithine (IS)", "2H3-Glutamate (IS)", "2H2-Citrulline (IS)", "2H4_13C5-Arginine (IS)",
                                  "13C6-Tyrosine (IS)")
internal_stand_pos_selection <- c("2H4-Alanine (IS)", "13C6-Phenylalanine (IS)", "2H4_13C5-Arginine (IS)", "2H3-Propionylcarnitine (IS)",
                                  "2H9-Isovalerylcarnitine (IS)")
internal_stand_sum_selection <- c("2H8-Valine (IS)", "2H3-Leucine (IS)", "2H3-Glutamate (IS)", "2H4_13C5-Arginine (IS)",
                                  "13C6-Tyrosine (IS)")

# add minimal intensity lines based on matrix (DBS or Plasma) and machine mode (neg, pos, sum)
if (dims_matrix == "DBS") {
  add_min_intens_lines <- TRUE
  hline_data_neg <-
    data.frame(
      int_line = c(15000, 200000, 130000, 18000, 50000),
      HMDB_name = internal_stand_neg_selection
    )
  hline_data_pos <-
    data.frame(
      int_line = c(150000, 3300000, 1750000, 150000, 270000),
      HMDB_name = internal_stand_pos_selection
    )
  hline_data_sum <-
    data.frame(
      int_line = c(1300000, 2500000, 500000, 1800000, 1400000),
      HMDB_name = internal_stand_sum_selection
    )
} else if (dims_matrix == "Plasma") {
  add_min_intens_lines <- TRUE
  hline_data_neg <-
    data.frame(
      int_line = c(6500, 100000, 75000, 7500, 25000),
      HMDB_name = internal_stand_neg_selection
    )
  hline_data_pos <-
    data.frame(
      int_line = c(85000, 1000000, 425000, 70000, 180000),
      HMDB_name = internal_stand_pos_selection
    )
  hline_data_sum <-
    data.frame(
      int_line = c(700000, 1250000, 150000, 425000, 300000),
      HMDB_name = internal_stand_sum_selection
    )
} else {
  add_min_intens_lines <- FALSE
}

# bar plots with a selection of IS
plot_width <- 8 + 0.2 * sample_count
plot_height <- plot_width / 2.5

is_neg_selection <- subset(internal_stand_neg, HMDB_name %in% internal_stand_neg_selection)
is_pos_selection <- subset(internal_stand_pos, HMDB_name %in% internal_stand_pos_selection)
is_sum_selection <- subset(internal_stand_summed, HMDB_name %in% internal_stand_sum_selection)

# bar plot either with or without minimal intensity lines
if (add_min_intens_lines) {
  save_internal_standard_plot(is_neg_selection, "barplot", "Interne Standaard (Neg)", outdir, "IS_bar_select_neg", plot_width, plot_height, hline_data_neg)
  save_internal_standard_plot(is_pos_selection, "barplot", "Interne Standaard (Pos)", outdir, "IS_bar_select_pos", plot_width, plot_height, hline_data_pos)
  save_internal_standard_plot(is_sum_selection, "barplot", "Interne Standaard (Sum)", outdir, "IS_bar_select_sum", plot_width, plot_height, hline_data_sum)
} else {
  save_internal_standard_plot(is_neg_selection, "barplot", "Interne Standaard (Neg)", outdir, "IS_bar_select_neg", plot_width, plot_height)
  save_internal_standard_plot(is_pos_selection, "barplot", "Interne Standaard (Pos)", outdir, "IS_bar_select_pos", plot_width, plot_height)
  save_internal_standard_plot(is_sum_selection, "barplot", "Interne Standaard (Sum)", outdir, "IS_bar_select_sum", plot_width, plot_height)
}

## line plots with a selection of IS
plot_width <- 8 + 0.2 * sample_count
plot_height <- plot_width / 2.0

save_internal_standard_plot(is_neg_selection, "lineplot", "Interne Standaard (Neg)", outdir, "IS_line_select_neg", plot_width, plot_height)
save_internal_standard_plot(is_pos_selection, "lineplot", "Interne Standaard (Pos)", outdir, "IS_line_select_pos", plot_width, plot_height)
save_internal_standard_plot(is_sum_selection, "lineplot", "Interne Standaard (Sum)", outdir, "IS_line_select_sum", plot_width, plot_height)

### POSITIVE CONTROLS CHECK
# these positive controls need to be in the samplesheet, in order to make the positive_control.RData file
# Positive control samples all have the format P1002.x, P1003.x and P1005.x (where x is a number)

column_list <- colnames(outlist)
patterns <- c("^(P1002\\.)[[:digit:]]+_", "^(P1003\\.)[[:digit:]]+_", "^(P1005\\.)[[:digit:]]+_")
positive_controls_index <- grepl(pattern = paste(patterns, collapse = "|"), column_list)
positive_control_list <- column_list[positive_controls_index]

if (z_score == 1) {
  # find if one or more positive control samples are missing
  pos_contr_warning <- c()
  if (all(sapply(c("^P1002", "^P1003", "^P1005"), 
                 function(x) any(grepl(x, positive_control_list))))) {
    cat("All three positive controls are present")
  } else {
    pos_contr_warning <- paste("positive controls list is not complete. Only",
                                paste(positive_control_list, collapse = ", "), "is/are present")
  }
  # you need all positive control samples, thus starting the script only if all are available
  if (length(positive_control_list) > 0) {
    # make positive control excel with specific HMDB_codes in combination with specific control samples
    positive_control <- NULL
    for (pos_ctrl in positive_control_list) {
      if (any(grepl("^P1002", pos_ctrl))) {
        pa_sample_name <- positive_control_list[grepl("P1002", positive_control_list)]
        pa_codes <- c("HMDB0000824", "HMDB0000725", "HMDB0000123")
        pa_names <- c("Propionylcarnitine", "Propionylglycine", "Glycine")
        pa_data <- get_pos_ctrl_data(outlist, pa_sample_name, pa_codes, pa_names)
        positive_control <- rbind(positive_control, pa_data)
      }
      if (any(grepl("^P1003", pos_ctrl))) {
        pku_sample_name <- positive_control_list[grepl("P1003", positive_control_list)]
        pku_codes <- c("HMDB0000159")
        pku_names <- c("L-Phenylalanine")
        pku_data <- get_pos_ctrl_data(outlist, pku_sample_name, pku_codes, pku_names)
        positive_control <- rbind(positive_control, pku_data)
      }
      if (any(grepl("^P1005", pos_ctrl))) {
        lpi_sample_name <- positive_control_list[grepl("P1005", positive_control_list)]
        lpi_codes <- c("HMDB0000904", "HMDB0000641", "HMDB0000182")
        lpi_names <- c("Citrulline", "L-Glutamine", "L-Lysine")
        lpi_data <- get_pos_ctrl_data(outlist, lpi_sample_name, lpi_codes, lpi_names)
        positive_control <- rbind(positive_control, lpi_data)
      }
    }

    positive_control$Zscore <- as.numeric(positive_control$Zscore)
    # extra information added to excel for future reference. made in beginning of this script
    positive_control$Matrix <- dims_matrix
    positive_control$Rundate <- rundate
    positive_control$Project <- project

    # Save results
    save(positive_control, file = paste0(outdir, "/", project, "_positive_control.RData"))
    # round the Z-scores to 2 digits
    positive_control$Zscore <- round_df(positive_control$Zscore, 2)
    write.xlsx(positive_control, file = paste0(outdir, "/", project, "_positive_control.xlsx"),
               sheetName = "Sheet1", colNames = TRUE, rowNames = TRUE, append = FALSE)
  }
  if (length(pos_contr_warning) != 0) {
    write.table(pos_contr_warning, file = paste(outdir, "positive_controls_warning.txt", sep = "/"),
                rowNames = FALSE, colNames = FALSE, quote = FALSE)
  }
}

### SST components output ####

# Internal standards lists, calculate coefficients of variation
if ("plots" %in% colnames(internal_stand_list)) {
  intensity_col_ids <- 2:(which(colnames(internal_stand_list) == "HMDB_name") - 1)
} else {
  intensity_col_ids <- 1:(which(colnames(internal_stand_list) == "HMDB_name") - 1)
}

internal_stand_list_intensities <- get_is_intensities(internal_stand_list, int_cols = intensity_col_ids)
internal_stand_neg_intensities <- get_is_intensities(outlist_tot_neg, is_codes = internal_stand_codes)
internal_stand_pos_intensities <- get_is_intensities(outlist_tot_pos, is_codes = internal_stand_codes)

# SST components.
sst_comp <- read.csv(sst_components_file, header = TRUE, sep = "\t")
sst_list <- outlist %>% filter(HMDB_code %in% sst_comp$HMDB_ID)
sst_colnrs <- grep("P1001", colnames(sst_list))

if (length(sst_colnrs) > 0) {
  sst_list_intensities <- sst_list[, sst_colnrs]
  control_col_ids <- grep(control_label, colnames(sst_list), fixed = TRUE)
  control_list_intensities <- sst_list[, control_col_ids]
  control_list_cv <- calculate_coefficient_of_variation(control_list_intensities)
  sst_list_intensities <- cbind(sst_list_intensities, CV_controls = control_list_cv[, "CV_perc"])
} else {
  sst_list_intensities <- sst_list[, intensity_col_ids]
}
for (col_nr in seq_len(ncol(sst_list_intensities))) {
  sst_list_intensities[, col_nr] <- as.numeric(sst_list_intensities[, col_nr])
  if (grepl("Zscore", colnames(sst_list_intensities)[col_nr])) {
    sst_list_intensities[, col_nr] <- round(sst_list_intensities[, col_nr], 2)
  } else {
    sst_list_intensities[, col_nr] <- round(sst_list_intensities[, col_nr])
  }
}
sst_list_intensities <- cbind(SST_comp_name = sst_list$HMDB_name, sst_list_intensities)

# Create Excel file
wb <- createWorkbook("IS_SST")
addWorksheet(wb, "Internal Standards")
openxlsx::writeData(wb, sheet = 1, internal_stand_list_intensities)
setColWidths(wb, 1, cols = 1, widths = 24)
addWorksheet(wb, "IS pos")
openxlsx::writeData(wb, sheet = 2, internal_stand_pos_intensities)
setColWidths(wb, 2, cols = 1, widths = 24)
addWorksheet(wb, "IS neg")
openxlsx::writeData(wb, sheet = 3, internal_stand_neg_intensities)
setColWidths(wb, 3, cols = 1, widths = 24)
addWorksheet(wb, "SST components")
openxlsx::writeData(wb, sheet = 4, sst_list_intensities)
setColWidths(wb, 4, cols = 1:3, widths = 24)
xlsx_name <- paste0(outdir, "/", project, "_IS_SST.xlsx")
openxlsx::saveWorkbook(wb, xlsx_name, overwrite = TRUE)
rm(wb)


### MISSING M/Z CHECK
# check the outlist_identified_(negative/positive).RData files for missing m/z values and save to file
# Load the outlist_identified files + remove the loaded files
load(paste0(outdir, "/outlist_identified_negative.RData"))
mzmed_pgrp_ident_neg <- outlist_ident$mzmed.pgrp
load(paste0(outdir, "/outlist_identified_positive.RData"))
mzmed_pgrp_ident_pos <- outlist_ident$mzmed.pgrp
rm(outlist_ident)

# Check for missing mz values, if present returned with vector of missing mz values
mz_missing_neg <- check_missing_mz(mzmed_pgrp_ident_neg, "Negative")
mz_missing_pos <- check_missing_mz(mzmed_pgrp_ident_pos, "Positive")

# Write both scanmodes to missing_mz_warning file
lapply(c(mz_missing_neg, mz_missing_pos), write, file = paste0(outdir, "/missing_mz_warning.txt"), append = TRUE, ncolumns = 1000)
