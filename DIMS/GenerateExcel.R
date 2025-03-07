## adapted from 13-excelExport.R

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
imagesize_multiplier <- 2

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
rlvnc_in_list <- rlvnc %>%
  filter(HMDB_key %in% rownames(outlist)[peaks_in_list]) %>%
  rename(sec_HMBD_ID_rlvnc = sec_HMDB_ID)
outlist <- cbind(outlist[peaks_in_list, ], as.data.frame(rlvnc_in_list))

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
wb_intensities <- openxlsx::createWorkbook("SinglePatient")
openxlsx::addWorksheet(wb_intensities, sheetname)

wb_helix_intensities <- openxlsx::createWorkbook("SinglePatient")
openxlsx::addWorksheet(wb_helix_intensities, sheetname)

# Add Z-scores and create plots
if (z_score == 1) {
  # add a column for plots
  outlist <- cbind(plots = NA, outlist)
  # two columns will be added for mean and stdev of controls; Z-scores start at ncol + 3
  startcol <- ncol(outlist) + 3

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

  # calculate Z-scores
  outlist <- calculate_zscores(outlist, "_Zscore", control_intensities, NULL, intensity_col_ids, startcol)

  # calculate robust Z-scores
  outlist_robustZ <- calculate_zscores(outlist, "_RobustZscore", control_col_idx, perc, intensity_col_ids, startcol)

  # calculate Z-scores after removal of outliers in Control samples with grubbs test
  outlist_nooutliers <- calculate_zscores(outlist, "_OutlierRemovedZscore", control_col_idx, outlier_threshold, intensity_col_ids, startcol)

  # output metabolites filtered on relevance
  save_to_rdata_and_txt(outlist, "AdductSums_filtered_Zscores")
  # output filtered metabolites with robust scaled Zscores
  save_to_rdata_and_txt(outlist_robustZ, "AdductSums_filtered_robustZ")
  # output filtered metabolites after removal of outliers
  save_to_rdata_and_txt(outlist_nooutliers, "AdductSums_filtered_outliersremovedZ")

  # get the IDs of the patients and sort
  patient_ids <- unique(gsub("\\.[0-9]*", "", patient_columns))
  patient_ids <- patient_ids[order(nchar(patient_ids), patient_ids)]

  # get Helix IDs for extra Excel file
  metabolite_files <- list.files(path = paste(path_metabolite_groups, "Diagnostics", sep = "/"),
                                 pattern = "*.txt", full.names = FALSE, recursive = FALSE)
  metab_df_helix <- NULL
  for (file_index in seq_along(metabolite_files)) {
    infile <- metabolite_files[file_index]
    metab_list <- read.table(paste(path_metabolite_groups, "Diagnostics", infile, sep = "/"),
                             sep = "\t", header = TRUE, quote = "")
    metab_df_helix <- rbind(metab_df_helix, metab_list)
  }
  # get Helix metabolites and unique HMDB IDs and remove ratio HMDBs containing A or L
  metab_df_helix <- metab_df_helix %>%
    filter(Helix == "ja") %>%
    select(c(HMDB_code, HMDB_name)) %>%
    rename(H_Name = HMDB_name)
  metab_list_helix <- unique(metab_df_helix$HMDB_code)
  metab_list_helix <- grep("[AL]", metab_list_helix, value = TRUE, invert = TRUE)

  outlist_helix <- outlist %>%
    filter(HMDB_key %in% metab_list_helix) %>%
    left_join(., metab_df_helix, by = join_by(HMDB_code == HMDB_code)) %>%
    select(-c(HMDB_key, sec_HMBD_ID_rlvnc, name, relevance, descr, origin, fluids, tissue, disease, pathway),
           -all_of(control_col_idx), -all_of(patient_col_idx)) %>%
    relocate(c(HMDB_code, H_Name, avg.ctrls, sd.ctrls), .after = plots) %>%
    relocate(c(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID), .after = last_col()) %>%
    rename(Name = H_Name)

  for (row_index in seq_len(nrow(outlist))) {
    # get HMDB ID
    hmdb_name <- rownames(outlist[row_index, ])

    # get intensities of controls and patient for a metabolite, get intensity columns,
    # pivot to long format, arrange Samples nummerically, change Sample names, get group size and
    # set Intensities to numeric.
    intensities <- outlist %>%
      slice(row_index) %>%
      select(all_of(intensity_col_ids)) %>%
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "Samples", values_to = "Intensities") %>%
      arrange(nchar(Samples)) %>%
      mutate(Samples = gsub("\\..*", "", Samples),
             Samples = gsub("(C).*", "\\1", Samples),
             Intensities = as.numeric(Intensities),
             type = ifelse(Samples == "C", "Control", "Patients")) %>%
      group_by(Samples) %>%
      mutate(group_size = n()) %>%
      ungroup()

    # set plot width to 40 times the number of samples
    plot_width <- length(unique(intensities$Samples)) * 40

    tmp_png <- "plot.png"
    png(filename = tmp_png, width = plot_width, height = 280)

    plot.new()
    # plot intensities for the controls and patients, use boxplot if group size is above 2, otherwise use a dash/line
    ggplot(intensities, aes(Samples, Intensities)) + geom_boxplot(data = subset(intensities, group_size > 2), aes(fill = type)) +
      theme_bw() +
      geom_point(data = subset(intensities, group_size <= 2), shape = "-", size = 5, aes(colour = type, fill = type)) +
      scale_fill_manual(values = c("green", "#930000")) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90), axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 10), axis.text = element_text(size = 7)) +
      ggtitle(hmdb_name)
    dev.off()

    # place the plot in the Excel file
    openxlsx::insertImage(wb_intensities,
                          sheetname,
                          tmp_png,
                          startRow = row_index + 1,
                          startCol = 1,
                          height = 560,
                          width = plot_width,
                          units = "px")

    if (hmdb_name %in% metab_list_helix) {
      openxlsx::insertImage(wb_helix_intensities,
                            sheetname,
                            tmp_png,
                            startRow = row_index + 1,
                            startCol = 1,
                            height = 560,
                            width = plot_width,
                            units = "px")
    }
  }
  wb_intensities <- set_row_height_col_width_wb(wb_intensities, sheetname, outlist, plot_width, plots_present = TRUE)

  wb_helix_intensities <- set_row_height_col_width_wb(wb_helix_intensities, sheetname, outlist_helix, plot_width, plots_present = TRUE)
  openxlsx::writeData(wb_helix_intensities, sheet = 1, outlist_helix, startCol = 1)
  openxlsx::saveWorkbook(wb_helix_intensities, paste0(outdir, "/Helix_", project, ".xlsx"), overwrite = TRUE)
  rm(wb_helix_intensities)
  # reorder outlist for Excel file
  outlist <- outlist %>% 
    relocate(c(HMDB_code, HMDB_name, avg.ctrls, sd.ctrls), .after = plots) %>%
    relocate(all_of(grep("_Zscore", colnames(outlist))), .after = sd.ctrls) %>% 
    relocate(all_of(c(colnames(control_intensities), patient_columns)), .after = last_col())
} else {
  wb_intensities <- set_row_height_col_width_wb(wb_intensities, sheetname, outlist, plot_width = NULL, plots_present = FALSE)
}

# write Excel file
openxlsx::writeData(wb_intensities, sheet = 1, outlist, startCol = 1)
openxlsx::saveWorkbook(wb_intensities, paste0(outdir, "/", project, ".xlsx"), overwrite = TRUE)
rm(wb_intensities)

# save outlist for GenerateQC step
save(outlist, file = "outlist.RData")
