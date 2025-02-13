## adapted from 13-excelExport.R

# load required packages
library("ggplot2")
library("reshape2")
library("openxlsx")
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

project <- cmd_args[1]
hmdb_rlvc_file <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])
export_scripts_dir <- cmd_args[4]

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

# create a directory for plots in project directory
dir.create(paste0(outdir, "/plots"), showWarnings = FALSE)
plot_dir <- paste0(outdir, "/plots/adducts")
dir.create(plot_dir, showWarnings = FALSE)

# load HMDB rlvnc table
load(hmdb_rlvc_file)

# load outlist object
load("AdductSums_combined.RData")

# Filter for biological relevance
peaks_in_list <- which(rownames(outlist) %in% rlvnc$HMDB_key)
rlvnc_in_list <- rlvnc %>% filter(HMDB_key %in% rownames(outlist)[peaks_in_list])
rlvnc_in_list <- rlvnc_in_list %>% rename(sec_HMBD_ID_rlvnc = sec_HMDB_ID)
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
filelist <- "AllPeakGroups"
wb_intensities <- createWorkbook("SinglePatient")
addWorksheet(wb_intensities, filelist)

# Add Z-scores and create plots
if (z_score == 1) {
  # add a column for plots
  outlist <- cbind(plots = NA, outlist)
  # two columns will be added for mean and stdev of controls; Z-scores start at ncol + 3
  startcol <- ncol(outlist) + 3

  # Get columns with control intensities
  control_intensity_cols <- get_intensities_cols(outlist, control_label)
  control_col_ids <- control_intensity_cols$col_ids
  control_columns <- control_intensity_cols$columns

  # Get columns with patient intensities
  patient_intensity_cols <- get_intensities_cols(outlist, case_label)
  patient_col_ids <- patient_intensity_cols$col_ids
  patient_columns <- patient_intensity_cols$columns

  intensity_col_ids <- c(control_col_ids, patient_col_ids)

  # if there are any intensities of 0 left, set them to NA for stats
  outlist[, intensity_col_ids][outlist[, intensity_col_ids] == 0] <- NA

  # calculate Z-scores
  outlist <- calculate_zscores(outlist, "_Zscore", control_columns, NULL, intensity_col_ids, startcol)

  # calculate robust Z-scores
  outlist_robustZ <- calculate_zscores(outlist, "_RobustZscore", control_col_ids, perc, intensity_col_ids, startcol)

  # calculate Z-scores after removal of outliers in Control samples with grubbs test
  outlist_nooutliers <- calculate_zscores(outlist, "_OutlierRemovedZscore", control_col_ids, outlier_threshold, intensity_col_ids, startcol)

  # output metabolites filtered on relevance
  save_to_rdata_and_txt(outlist, "AdductSums_filtered_Zscores")
  # output filtered metabolites with robust scaled Zscores
  save_to_rdata_and_txt(outlist_robustZ, "AdductSums_filtered_robustZ")
  # output filtered metabolites after removal of outliers
  save_to_rdata_and_txt(outlist_nooutliers, "AdductSums_filtered_outliersremovedZ")

  # get the IDs of the patients and sort
  patient_ids <- unique(sapply(strsplit(colnames(patient_columns), ".", fixed = TRUE), "[", 1))
  patient_ids <- patient_ids[order(nchar(patient_ids), patient_ids)]

  # for every row make a plot and insert into excel
  temp_png <- NULL
  for (p in seq_len(nrow(outlist))) {
    # get HMDB ID
    hmdb_name <- rownames(outlist[p, ])

    # get intensities of controls and patient for a metabolite, get intensity columns,
    # pivot to long format, arrange Samples nummerically, change Sample names, get group size and
    # set Intensities to numeric.
    intensities <- outlist %>%
      slice(p) %>%
      select(intensity_col_ids) %>%
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "Samples", values_to = "Intensities") %>%
      arrange(nchar(Samples)) %>%
      group_by(Samples) %>%
      mutate(Samples = gsub("\\..*", "", Samples),
             Samples = gsub("(C).*", "\\1", Samples),
             group_size = n()) %>%
      ungroup() %>%
      as.numeric(Intensities)

    plot.new()

    # plot intensities for the controls and patients, use boxplot if group size is above 2, otherwise use a dash/line
    ggplot(intensities, aes(Samples, Intensities)) + geom_boxplot(data = subset(intensities, group_size > 2)) +
      theme_bw() +
      geom_point(data = subset(intensities, group_size <= 2), shape = "-", size = 5, aes(colour = type, fill = type)) +
      scale_fill_manual(values = c("green", "#930000")) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90),
            axis.title = element_blank(), plot.title = element_text(hjust = 0.5)) +
      ggtitle(hmdb_name)

    # set plot width to 40 times the number of samples
    plot_width <- length(unique(intensities$Samples)) * 40

    # place the plot in the Excel file
    openxlsx::insertImage(wb_intensities,
                          filelist,
                          file_png,
                          startRow = p + 1,
                          startCol = 1,
                          height = 560,
                          width = plot_width,
                          units = "px")
  }

  openxlsx::setRowHeights(wb_intensities, filelist, rows = c(seq_len(1, nrow(outlist) + 1)), heights = cell_dim[1] / 4)
  openxlsx::setColWidths(wb_intensities, filelist, cols = c(seq_len(2, ncol(outlist))), widths = 20)
} else {
  openxlsx::setRowHeights(wb_intensities, filelist, rows = c(seq_len(1, nrow(outlist))), heights = 18)
  openxlsx::setColWidths(wb_intensities, filelist, cols = c(seq_len(1, ncol(outlist))), widths = 20)
}

# write Excel file
openxlsx::writeData(wb_intensities, sheet = 1, outlist, startCol = 1)
xlsx_name <- paste0(outdir, "/", project, ".xlsx")
openxlsx::saveWorkbook(wb_intensities, xlsx_name, overwrite = TRUE)
rm(wb_intensities)

# save outlist for GenerateQC step
save(outlist, file = "outlist.RData")
