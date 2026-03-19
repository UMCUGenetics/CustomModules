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
sst_components_file <- cmd_args[4]
export_scripts_dir <- cmd_args[5]

outdir <- "./"

# load in function scripts
source(paste0(export_scripts_dir, "generate_qc_output_functions.R"))

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

#### CHECK NUMBER OF CONTROLS ####
file_name <- "Check_number_of_controls.txt"
min_num_controls <- 25
check_number_of_controls(outlist, min_num_controls, file_name)

#### INTERNAL STANDARDS ####
is_list <- outlist[grep("Internal standard", outlist[, "relevance"], fixed = TRUE), ]
is_codes <- rownames(is_list)

# check if there is data present for all the samples that the pipeline started with
sample_names_nodata <- setdiff(names(repl_pattern), names(is_list))
if (length(sample_names_nodata) == 0) {
  sample_names_nodata <- "none"
}
write.table(sample_names_nodata,
  file = paste(outdir, "sample_names_nodata.txt", sep = "/"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)
if (!is.null(sample_names_nodata)) {
  for (sample_name in sample_names_nodata) {
    repl_pattern[[sample_name]] <- NULL
  }
}

# Only continue with patients (columns) that are in both pos and neg, so patients that are in both
samples_both_modes <- intersect(colnames(outlist_tot_neg), colnames(outlist_tot_pos))
outlist_tot_neg <- outlist_tot_neg[, samples_both_modes]
outlist_tot_pos <- outlist_tot_pos[, samples_both_modes]

# Retrieve IS summed adducts
is_summed <- get_internal_standards(is_list, "summed", repl_pattern, dims_matrix, rundate, project)
# Retrieve IS positive mode
is_pos <- get_internal_standards(is_list, "pos", outlist_tot_pos, dims_matrix, rundate, project)
# Retrieve IS negative mode
is_neg <- get_internal_standards(is_list, "neg", outlist_tot_neg, dims_matrix, rundate, project)

# Save results
save(is_pos, is_neg, is_summed, file = paste0(outdir, "/", project, "_IS_results.RData"))

# number of samples, for plotting length and width
sample_count <- length(repl_pattern)

# change the order of the x-axis summed plots to a natural sorted one
sample_naturalorder <- unique(as.character(is_summed$Sample))
sample_naturalorder <- stringr::str_sort(sample_naturalorder, numeric = TRUE)
is_summed$Sample_level <- factor(is_summed$Sample, levels = c(sample_naturalorder))
is_pos$Sample_level <- factor(is_pos$Sample, levels = c(sample_naturalorder))
is_neg$Sample_level <- factor(is_neg$Sample, levels = c(sample_naturalorder))

## bar plots with all IS
plot_width <- 9 + 0.35 * sample_count
plot_height <- plot_width / 2.5

save_internal_standard_plot(
  is_neg, "barplot", "Interne Standaard (Neg)", outdir,
  "IS_bar_all_neg", plot_width, plot_height
)
save_internal_standard_plot(
  is_pos, "barplot", "Interne Standaard (Pos)", outdir,
  "IS_bar_all_pos", plot_width, plot_height
)
save_internal_standard_plot(
  is_summed, "barplot", "Interne Standaard (Summed)", outdir,
  "IS_bar_all_sum", plot_width, plot_height
)

## Line plots with all IS
plot_width <- 8 + 0.2 * sample_count
plot_height <- plot_width / 2.5

save_internal_standard_plot(
  is_neg, "lineplot", "Interne Standaard (Neg)",
  outdir, "IS_line_all_neg", plot_width, plot_height
)
save_internal_standard_plot(
  is_pos, "lineplot", "Interne Standaard (Pos)",
  outdir, "IS_line_all_pos", plot_width, plot_height
)
save_internal_standard_plot(
  is_summed, "lineplot", "Interne Standaard (Sum)",
  outdir, "IS_line_all_sum", plot_width, plot_height
)

## bar plots with a selection of IS
is_neg_selection <- c(
  "2H2-Ornithine (IS)", "2H3-Glutamate (IS)",
  "2H2-Citrulline (IS)", "2H4_13C5-Arginine (IS)",
  "13C6-Tyrosine (IS)"
)
is_pos_selection <- c(
  "2H4-Alanine (IS)", "13C6-Phenylalanine (IS)",
  "2H4_13C5-Arginine (IS)", "2H3-Propionylcarnitine (IS)",
  "2H9-Isovalerylcarnitine (IS)"
)
is_sum_selection <- c(
  "2H8-Valine (IS)", "2H3-Leucine (IS)",
  "2H3-Glutamate (IS)", "2H4_13C5-Arginine (IS)",
  "13C6-Tyrosine (IS)"
)

# define threshold for acceptance of selected internal standards
threshold_is_dbs_neg <- c(15000, 200000, 130000, 18000, 50000)
threshold_is_dbs_pos <- c(150000, 3300000, 1750000, 150000, 270000)
threshold_is_dbs_sum <- c(1300000, 2500000, 500000, 1800000, 1400000)
threshold_is_pl_neg  <- c(70000, 700000, 700000, 65000, 350000)
threshold_is_pl_pos  <- c(1500000, 9000000, 3000000, 400000, 700000)
threshold_is_pl_sum  <- c(8000000, 12500000, 2500000, 3000000, 4000000)

# add minimal intensity lines based on matrix (DBS or Plasma) and machine mode (neg, pos, sum)
if (dims_matrix == "DBS") {
  add_min_intens_lines <- TRUE
  hline_data_neg <-
    data.frame(
      int_line = threshold_is_dbs_neg,
      HMDB_name = is_neg_selection
    )
  hline_data_pos <-
    data.frame(
      int_line = threshold_is_dbs_pos,
      HMDB_name = is_pos_selection
    )
  hline_data_sum <-
    data.frame(
      int_line = threshold_is_dbs_sum,
      HMDB_name = is_sum_selection
    )
} else if (dims_matrix == "Plasma") {
  add_min_intens_lines <- TRUE
  hline_data_neg <-
    data.frame(
      int_line = threshold_is_pl_neg,
      HMDB_name = is_neg_selection
    )
  hline_data_pos <-
    data.frame(
      int_line = threshold_is_pl_pos,
      HMDB_name = is_pos_selection
    )
  hline_data_sum <-
    data.frame(
      int_line = threshold_is_pl_sum,
      HMDB_name = is_sum_selection
    )
} else {
  add_min_intens_lines <- FALSE
}

# bar plots with a selection of IS
plot_width <- 8 + 0.2 * sample_count
plot_height <- plot_width / 2.5

is_neg_selection_subset <- subset(is_neg, HMDB_name %in% is_neg_selection)
is_pos_selection_subset <- subset(is_pos, HMDB_name %in% is_pos_selection)
is_sum_selection_subset <- subset(is_summed, HMDB_name %in% is_sum_selection)

# export txt file with samples with internal standard level below threshold
if (dims_matrix == "Plasma") {
  is_below_threshold_neg <- find_is_below_threshold(is_neg_selection_subset, threshold_is_pl_neg, is_neg_selection, "neg")
  is_below_threshold_pos <- find_is_below_threshold(is_pos_selection_subset, threshold_is_pl_pos, is_pos_selection, "pos")
  is_below_threshold_sum <- find_is_below_threshold(is_sum_selection_subset, threshold_is_pl_sum, is_sum_selection, "sum")
  is_below_threshold <- rbind(is_below_threshold_pos, is_below_threshold_neg, is_below_threshold_sum)
} else if (dims_matrix == "DBS") {
  is_below_threshold_neg <- find_is_below_threshold(is_neg_selection_subset, threshold_is_dbs_neg, is_neg_selection, "neg")
  is_below_threshold_pos <- find_is_below_threshold(is_pos_selection_subset, threshold_is_dbs_pos, is_pos_selection, "pos")
  is_below_threshold_sum <- find_is_below_threshold(is_sum_selection_subset, threshold_is_dbs_sum, is_neg_selection, "sum")
  is_below_threshold <- rbind(is_below_threshold_pos, is_below_threshold_neg, is_below_threshold_sum)
} else {
  # generate empty table
  is_below_threshold <- is_neg_selection_subset[0, ]
}

if (nrow(is_below_threshold) > 0) {
  write.table(is_below_threshold, 
	      file = "internal_standards_below_threshold.txt", 
	      row.names = FALSE, sep = "\t")
} else { 
  write.table("no internal standards are below threshold",
	      file = "internal_standards_below_threshold.txt",
	      row.names = FALSE, col.names = FALSE
	      )
}

# bar plot either with or without minimal intensity lines
if (add_min_intens_lines) {
  save_internal_standard_plot(
    is_neg_selection_subset, "barplot", "Interne Standaard (Neg)", outdir,
    "IS_bar_select_neg", plot_width, plot_height, hline_data_neg
  )
  save_internal_standard_plot(
    is_pos_selection_subset, "barplot", "Interne Standaard (Pos)", outdir,
    "IS_bar_select_pos", plot_width, plot_height, hline_data_pos
  )
  save_internal_standard_plot(
    is_sum_selection_subset, "barplot", "Interne Standaard (Sum)", outdir,
    "IS_bar_select_sum", plot_width, plot_height, hline_data_sum
  )
} else {
  save_internal_standard_plot(
    is_neg_selection_subset, "barplot", "Interne Standaard (Neg)", outdir,
    "IS_bar_select_neg", plot_width, plot_height
  )
  save_internal_standard_plot(
    is_pos_selection_subset, "barplot", "Interne Standaard (Pos)", outdir,
    "IS_bar_select_pos", plot_width, plot_height
  )
  save_internal_standard_plot(
    is_sum_selection_subset, "barplot", "Interne Standaard (Sum)", outdir,
    "IS_bar_select_sum", plot_width, plot_height
  )
}

## line plots with a selection of IS
plot_width <- 8 + 0.2 * sample_count
plot_height <- plot_width / 2.0

save_internal_standard_plot(
  is_neg_selection_subset, "lineplot", "Interne Standaard (Neg)", outdir,
  "IS_line_select_neg", plot_width, plot_height
)
save_internal_standard_plot(
  is_pos_selection_subset, "lineplot", "Interne Standaard (Pos)", outdir,
  "IS_line_select_pos", plot_width, plot_height
)
save_internal_standard_plot(
  is_sum_selection_subset, "lineplot", "Interne Standaard (Sum)", outdir,
  "IS_line_select_sum", plot_width, plot_height
)

### POSITIVE CONTROLS CHECK
# these positive controls need to be in the samplesheet, in order to make the positive_control.RData file
# Positive control samples all have the format P1002.x, P1003.x and P1005.x (where x is a number)

column_list <- colnames(outlist)
patterns <- c("^(P1002\\.)[[:digit:]]+_", "^(P1003\\.)[[:digit:]]+_", "^(P1005\\.)[[:digit:]]+_")
positive_controls_index <- grepl(pattern = paste(patterns, collapse = "|"), column_list)
positive_control_list <- column_list[positive_controls_index]

if (sum(positive_controls_index) > 0) {
  # find if one or more positive control samples are missing
  pos_contr_warning <- c()
  if (all(sapply(c("^P1002", "^P1003", "^P1005"),
                 function(x) any(grepl(x, positive_control_list))))) {
    pos_contr_warning <- "All three positive controls are present"
  } else {
    pos_contr_warning <- paste(
      "positive controls list is not complete. Only",
      paste(positive_control_list, collapse = ", "), " present"
    )
  }
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
    write.xlsx(positive_control,
      file = paste0(outdir, "/", project, "_positive_control.xlsx"),
      sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE
    )
  } 
  if (length(pos_contr_warning) == 0) {
    pos_contr_warning <- "No positive controls found"
  }
  write.table(pos_contr_warning,
    file = paste(outdir, "positive_controls_warning.txt", sep = "/"),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )
}

### SST components output ####

# Internal standards lists, calculate coefficients of variation
if ("plots" %in% colnames(is_list)) {
  intensity_col_ids <- 2:(which(colnames(is_list) == "HMDB_name") - 1)
} else {
  intensity_col_ids <- 1:(which(colnames(is_list) == "HMDB_name") - 1)
}

is_list_intensities <- get_is_intensities(is_list, int_cols = intensity_col_ids)
is_neg_intensities <- get_is_intensities(outlist_tot_neg, is_codes = is_codes)
is_pos_intensities <- get_is_intensities(outlist_tot_pos, is_codes = is_codes)

# SST components
sst_components <- read.csv(sst_components_file, header = TRUE, sep = "\t")
sst_metabolites_df <- outlist %>% filter(HMDB_code %in% sst_components$HMDB_ID)
sst_sample_column_index <- grep("P1001", colnames(sst_metabolites_df))

# Check if SST mix sample(s) are present
if (length(sst_sample_column_index) > 0) {
  # Get the SST intensities of the controls, calculate the coefficient of variation
  # and add to SST mix intensities
  sst_sample_intensities_df <- sst_metabolites_df[, sst_sample_column_index]
  control_col_ids <- grep(control_label, colnames(sst_metabolites_df), fixed = TRUE)
  control_sst_intensities_df <- sst_metabolites_df[, control_col_ids]
  control_sst_metabolites_cv <- calc_coefficient_of_variation(control_sst_intensities_df)
  sst_intensities_df <- cbind(sst_sample_intensities_df, CV_controls = control_sst_metabolites_cv[, "CV_perc"])
} else {
  # Use intensities when there is not SST mix sample added
  sst_intensities_df <- sst_metabolites_df[, intensity_col_ids]
}

sst_intensities_df <- as.data.frame(sst_intensities_df)
for (col_nr in seq_len(ncol(sst_intensities_df))) {
  # Change column type to numeric
  sst_intensities_df[, col_nr] <- as.numeric(sst_intensities_df[, col_nr])
  if (grepl("Zscore", colnames(sst_intensities_df)[col_nr])) {
    # Round numeric value of Z-score columns to 2 decimal places
    sst_intensities_df[, col_nr] <- round(sst_intensities_df[, col_nr], 2)
  } else {
    # Round numeric value of intensity columns to an intiger
    sst_intensities_df[, col_nr] <- round(sst_intensities_df[, col_nr])
  }
}
sst_intensities_df <- cbind(SST_comp_name = sst_metabolites_df$HMDB_name, sst_intensities_df)

# Create Excel file
wb <- createWorkbook("IS_SST")
addWorksheet(wb, "Internal Standards")
openxlsx::writeData(wb, sheet = 1, is_list_intensities)
setColWidths(wb, 1, cols = 1, widths = 24)
addWorksheet(wb, "IS pos")
openxlsx::writeData(wb, sheet = 2, is_pos_intensities)
setColWidths(wb, 2, cols = 1, widths = 24)
addWorksheet(wb, "IS neg")
openxlsx::writeData(wb, sheet = 3, is_neg_intensities)
setColWidths(wb, 3, cols = 1, widths = 24)
addWorksheet(wb, "SST components")
openxlsx::writeData(wb, sheet = 4, sst_intensities_df)
setColWidths(wb, 4, cols = 1:3, widths = 24)
xlsx_name <- paste0(outdir, "/", project, "_IS_SST.xlsx")
openxlsx::saveWorkbook(wb, xlsx_name, overwrite = TRUE)
rm(wb)

# generate text file for workflow completed mail for components with Z-score < 2
if (sum(grepl("P1001", colnames(sst_intensities_df))) > 0) {
  zscore_column <- grep("_Zscore", colnames(sst_intensities_df))[1]
  sst_intensities_df_qc <- sst_intensities_df[sst_intensities_df[, zscore_column] < 2, ]
  sst_intensities_df_qc <- select(sst_intensities_df_qc, -c("CV_controls"))
  write.table(sst_intensities_df_qc, file = paste(outdir, "sst_qc.txt", sep = "/"), row.names = FALSE, sep = "\t")
  # in case of an empty table, the column header doesn't need to appear in the mail
  if (nrow(sst_intensities_df_qc) == 0) {
    write.table("none", file = paste(outdir, "sst_qc.txt", sep = "/"), row.names = FALSE, col.names = FALSE)
  }
} else {
  write.table("no SST sample present", file = paste(outdir, "sst_qc.txt", sep = "/"), row.names = FALSE, col.names = FALSE)
}

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
lapply(c(mz_missing_neg, mz_missing_pos), write,
  file = paste0(outdir, "/missing_mz_warning.txt"),
  append = TRUE, ncolumns = 1000
)
