# load required packages
library("ggplot2")
library("reshape2")
library("openxlsx")
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

project <- cmd_args[1]
hmdb_rlvc_file <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])
export_scripts_dir <- cmd_args[4]
path_metabolite_groups <- cmd_args[5]

# load in function scripts
source(paste0(export_scripts_dir, "generate_excel_functions.R"))

# set the number of digits for floats
options(digits = 16)

# Initialise
plot <- TRUE
export <- TRUE
control_label <- "C"
case_label <- "P"

# setting outdir to export files to the working directory
outdir <- "./"
# percentage of outliers to remove from calculation of robust scaler
perc <- 5
# Z-score for removing outliers with grubbs test
outlier_threshold <- 2

# load HMDB rlvnc table
load(hmdb_rlvc_file)

# load outlist object
load("AdductSums_combined.RData")

# Filter for biological relevance
peaks_in_list <- which(rownames(outlist) %in% rlvnc$HMDB_key)
outlist_subset <- outlist[peaks_in_list, ]
outlist_subset$HMDB_key <- rownames(outlist_subset)
outlist <- outlist_subset %>%
  left_join(rlvnc %>% rename(sec_HMDB_ID_rlvnc = sec_HMDB_ID), by = "HMDB_key")
rownames(outlist) <- outlist$HMDB_key

# filter out all irrelevant HMDBs
outlist <- outlist %>%
  tibble::rownames_to_column("rowname") %>%
  filter(grepl("relevant|Onbekend|Internal", relevance)) %>%
  tibble::column_to_rownames("rowname")

# Add HMDB_code column with all the HMDB ID and sort on it
outlist <- cbind(outlist, "HMDB_code" = rownames(outlist))
outlist <- outlist[order(outlist[, "HMDB_code"]), ]

# Create excel
sheetname <- "AllPeakGroups"
wb_intensities_zscores <- openxlsx::createWorkbook("SinglePatient")
openxlsx::addWorksheet(wb_intensities_zscores, sheetname)

# Add Z-scores and create plots
if (z_score == 1) {
  dir.create(paste0(outdir, "/plots"), showWarnings = FALSE)
  wb_helix_zscores <- openxlsx::createWorkbook("SinglePatient")
  openxlsx::addWorksheet(wb_helix_zscores, sheetname)
  row_helix <- 2 # start on row 2 because of header
  # add a column for plots
  outlist <- cbind(plots = NA, outlist)
  # three columns will be added for mean, stdev and number of controls; Z-scores start at ncol + 4
  startcol <- ncol(outlist) + 4

  # Get columns with control intensities
  control_intensity_cols <- get_intensities_cols(outlist, control_label)
  control_col_idx <- control_intensity_cols$col_idx
  control_intensities <- control_intensity_cols$df_intensities

  # Get columns with patient intensities
  patient_intensity_cols <- get_intensities_cols(outlist, case_label)
  patient_col_idx <- patient_intensity_cols$col_idx
  patient_columns <- colnames(patient_intensity_cols$df_intensities)

  intensity_col_ids <- c(control_col_idx, patient_col_idx)

  # if there are any intensities of 0 left, set them to NA for stats
  outlist[, intensity_col_ids][outlist[, intensity_col_ids] == 0] <- NA

  # calculate robust Z-scores
  outlist_robust_zscore <- calculate_zscores(outlist, "_RobustZscore", control_col_idx, perc, intensity_col_ids, startcol)

  # calculate Z-scores after removal of outliers in Control samples with grubbs test
  outlist_nooutliers <- calculate_zscores(
    outlist, "_OutlierRemovedZscore", control_col_idx, outlier_threshold,
    intensity_col_ids, startcol
  )

  # calculate Z-scores
  outlist <- calculate_zscores(outlist, "_Zscore", control_intensities, NULL, intensity_col_ids, startcol)

  # output metabolites filtered on relevance
  save_to_rdata_and_txt(outlist, "AdductSums_filtered_Zscores")
  # output filtered metabolites with robust scaled Zscores
  save_to_rdata_and_txt(outlist_robust_zscore, "AdductSums_filtered_robustZ")
  # output filtered metabolites after removal of outliers
  save_to_rdata_and_txt(outlist_nooutliers, "AdductSums_filtered_outliersremovedZ")

  # use outlier-removed outlist for generating Excel file
  outlist <- outlist_nooutliers
  colnames(outlist) <- gsub("_OutlierRemovedZscore", "_Zscore", colnames(outlist))

  # save outlist for GenerateQC step
  save(outlist, file = "outlist.RData")

  # get Helix IDs for extra Excel file
  metabolite_files <- list.files(
    path = paste(path_metabolite_groups, "Diagnostics", sep = "/"),
    pattern = "*.txt", full.names = FALSE, recursive = FALSE
  )
  metab_df_helix <- NULL
  for (file_index in seq_along(metabolite_files)) {
    infile <- metabolite_files[file_index]
    metab_list <- read.table(paste(path_metabolite_groups, "Diagnostics", infile, sep = "/"),
      sep = "\t", header = TRUE, quote = ""
    )
    metab_df_helix <- rbind(metab_df_helix, metab_list)
  }
  # get Helix metabolites and unique HMDB IDs and remove ratio HMDBs containing A or L
  metab_df_helix <- metab_df_helix %>%
    filter(Helix == "ja") %>%
    select(c(HMDB_code, HMDB_name)) %>%
    rename(helix_name = HMDB_name)
  metab_list_helix <- unique(metab_df_helix$HMDB_code)
  metab_list_helix <- grep("[AL]", metab_list_helix, value = TRUE, invert = TRUE)

  outlist_helix <- outlist %>%
    filter(HMDB_key %in% metab_list_helix) %>%
    left_join(., metab_df_helix, by = join_by(HMDB_code == HMDB_code)) %>%
    select(
      -c(HMDB_key, sec_HMDB_ID_rlvnc, theormz_HMDB, name, relevance, descr, origin, fluids, tissue, disease, pathway),
      -all_of(control_col_idx), -all_of(patient_col_idx)
    ) %>%
    relocate(c(HMDB_code, helix_name, avg_ctrls, sd_ctrls), .after = plots) %>%
    relocate(c(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID), .after = last_col()) %>%
    rename(Name = helix_name)

  # Get intensity columns for controls and patients
  intensities_df <- outlist %>% select(HMDB_key, matches("^C|^P[0-9]"), -ends_with("_Zscore"))

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

    if (hmdb_id %in% metab_list_helix) {
      # Make separate plot for Helix Excel containing all samples

      start_row_index <- row_index + 1
      save_plot_to_excel_workbook(
        wb_helix_zscores,
        sheetname,
        intensities_df_long,
        "plots/plot_helix_",
        hmdb_id,
        plot_width,
        col_width,
        row_helix
      )
      row_helix <- row_helix + 1
    }

    # Remove postive controls and SST mix samples, (e.g. P1001, P1002, P1003, P1005)
    intensities_df_long <- intensities_df_long %>% filter(!grepl("^P[0-9]{4}$", Samples))

    start_row_index <- row_index + 1
    save_plot_to_excel_workbook(
      wb_intensities_zscores,
      sheetname,
      intensities_df_long,
      "plots/plot_",
      hmdb_id,
      plot_width,
      col_width,
      start_row_index
    )
  }
  wb_intensities <- set_row_height_col_width_wb(
    wb_intensities_zscores,
    sheetname,
    nrow(outlist),
    ncol(outlist),
    col_width,
    plots_present = TRUE
  )

  wb_helix_intensities <- set_row_height_col_width_wb(
    wb_helix_zscores,
    sheetname,
    nrow(outlist_helix),
    ncol(outlist_helix),
    col_width,
    plots_present = TRUE
  )
  openxlsx::writeData(wb_helix_intensities, sheet = 1, outlist_helix, startCol = 1)
  openxlsx::saveWorkbook(wb_helix_intensities, paste0(outdir, "/Helix_", project, ".xlsx"), overwrite = TRUE)
  rm(wb_helix_intensities)

  # reorder outlist for Excel file
  outlist <- outlist %>%
    relocate(c(HMDB_code, HMDB_name_all, theormz_HMDB, descr, avg_ctrls, sd_ctrls), .after = plots) %>%
    relocate(all_of(grep("_Zscore", colnames(outlist))), .after = sd_ctrls) %>%
    relocate(all_of(c(colnames(control_intensities), patient_columns)), .after = last_col())
} else {
  save(outlist, file = "outlist.RData")
  wb_intensities <- set_row_height_col_width_wb(
    wb_intensities_zscores,
    sheetname,
    nrow(outlist),
    ncol(outlist),
    plot_width = NULL,
    plots_present = FALSE
  )
  outlist <- outlist %>%
    relocate(c(HMDB_name, HMDB_name_all, HMDB_code, HMDB_ID_all, theormz_HMDB))
}

# write Excel file
openxlsx::writeData(wb_intensities_zscores, sheet = 1, outlist, startCol = 1)
openxlsx::saveWorkbook(wb_intensities_zscores, paste0(outdir, "/", project, ".xlsx"), overwrite = TRUE)
rm(wb_intensities_zscores)
unlink("plots", recursive = TRUE)
