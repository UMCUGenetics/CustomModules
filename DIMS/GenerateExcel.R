## adapted from 13-excelExport.R

# load required packages
library("ggplot2")
library("reshape2")
library("openxlsx")
library("loder")
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

init_file <- cmd_args[1]
project <- cmd_args[2]
dims_matrix <- cmd_args[3]
hmdb_file <- cmd_args[4]
z_score <- as.numeric(cmd_args[5])

round_df <- function(df, digits) {
  #' Round numbers to a set number of digits for numeric values
  #'
  #' @param df: Dataframe containing numeric values
  #' @param digits: Number of digits to round off to (integer)
  #'
  #' @return df: Dataframe with rounded numbers
  numeric_columns <- sapply(df, mode) == "numeric"
  df[numeric_columns] <- round(df[numeric_columns], digits)
  return(df)
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

# load information on samples
load(init_file)
# load the HMDB file with info on biological relevance of metabolites
load(hmdb_file)

# get current date
rundate <- Sys.Date()

# create a directory for plots in project directory
dir.create(paste0(outdir, "/plots"), showWarnings = FALSE)
plot_dir <- paste0(outdir, "/plots/adducts")
dir.create(plot_dir, showWarnings = FALSE)

# set the number of digits for floats
options(digits = 16)

# load positive and negative adduct sums
load("AdductSums_negative.RData")
outlist_neg_adducts_hmdb <- outlist.tot
load("AdductSums_positive.RData")
outlist_pos_adducts_hmdb <- outlist.tot
rm(outlist.tot)

# Only continue with patients (columns) that are in both pos and neg, so patients that are in both
tmp <- intersect(colnames(outlist_neg_adducts_hmdb), colnames(outlist_pos_adducts_hmdb))
outlist_neg_adducts_hmdb <- outlist_neg_adducts_hmdb[, tmp]
outlist_pos_adducts_hmdb <- outlist_pos_adducts_hmdb[, tmp]

# Find indexes of neg hmdb code that are also found in pos and vice versa
index_neg <- which(rownames(outlist_neg_adducts_hmdb) %in% rownames(outlist_pos_adducts_hmdb))
index_pos <- which(rownames(outlist_pos_adducts_hmdb) %in% rownames(outlist_neg_adducts_hmdb))

# Only continue with HMDB codes (rows) that were found in both positive mode and remove last column (hmdb_name)
tmp_pos <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ] %>% select(-c(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID))
tmp_pos_info <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ]  %>% select(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID)
# tmp_hmdb_name_pos <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], dim(outlist_pos_adducts_hmdb)[2]]
tmp_pos_left <- outlist_pos_adducts_hmdb[-index_pos, ]
# same for negative mode
tmp_neg <- outlist_neg_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], ] %>% select(-c(HMDB_name, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID))
tmp_neg_left <- outlist_neg_adducts_hmdb[-index_neg, ]
print("combine pos + neg")
# Combine positive and negative numbers and paste back HMDB column
tmp <- apply(tmp_pos, 2, as.numeric) + apply(tmp_neg, 2, as.numeric)
rownames(tmp) <- rownames(tmp_pos)
tmp <- cbind(tmp, tmp_pos_info)
outlist <- rbind(tmp, tmp_pos_left, tmp_neg_left)
outlist <- outlist %>% arrange(rownames(outlist))
print("Add biological relevance")

# Filter for biological relevance
# peaks_in_list <- which(rownames(outlist) %in% rownames(rlvnc))
peaks_in_list <- which(rownames(outlist) %in% rlvnc$HMDB_key)
rlvnc_in_list <- rlvnc %>% filter(HMDB_key %in% rownames(outlist)[peaks_in_list])
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
wb <- createWorkbook("SinglePatient")
addWorksheet(wb, filelist)

# Add Z-scores and create plots
if (z_score == 1) {
  # add a column for plots
  outlist <- cbind(plots = NA, outlist)
  # two columns will be added for mean and stdev of controls; Z-scores start at ncol + 3
  startcol <- ncol(outlist) + 3

  # Get columns with control intensities
  control_col_ids <- grep(control_label, colnames(outlist), fixed = TRUE)
  control_columns <- as.data.frame(outlist[, control_col_ids])
  colnames(control_columns) <- colnames(outlist)[control_col_ids]

  # Get columns with patient intensities
  patient_col_ids <- grep(case_label, colnames(outlist), fixed = TRUE)
  patient_columns <- as.data.frame(outlist[, patient_col_ids])
  colnames(patient_columns) <- colnames(outlist)[patient_col_ids]

  intensity_col_ids <- c(control_col_ids, patient_col_ids)

  # if there are any intensities of 0 left, set them to NA for stats
  outlist[, intensity_col_ids][outlist[, intensity_col_ids] == 0] <- NA

  # save outlist as it is and use it to calculate robust scaler
  outlist_noZ <- outlist

  # calculate mean and sd for Control group
  outlist$avg.ctrls <- apply(control_columns, 1, function(x) mean(as.numeric(x), na.rm = TRUE))
  outlist$sd.ctrls  <- apply(control_columns, 1, function(x) sd(as.numeric(x), na.rm = TRUE))

  # Make and add columns with zscores
  colnames_z <- NULL
  for (i in intensity_col_ids) {
    cname <- colnames(outlist)[i]
    colnames_z <- c(colnames_z, paste(cname, "Zscore", sep = "_"))
    zscores_1col <- (as.numeric(as.vector(unlist(outlist[, i]))) - outlist$avg.ctrls) / outlist$sd.ctrls
    outlist <- cbind(outlist, zscores_1col)
  }
  colnames(outlist)[startcol:ncol(outlist)] <- colnames_z

  # calculate robust scaler (Zscores minus outliers in Controls)
  outlist_noZ$avg.ctrls <- 0
  outlist_noZ$sd.ctrls  <- 0

  # only calculate robust Z-scores if there are enough Controls
  if (length(control_col_ids) > 10) {
    for (metabolite_index in 1:nrow(outlist)) {
      outlist_noZ$avg.ctrls[metabolite_index] <- mean(robust_scaler(outlist_noZ[metabolite_index, control_col_ids],
                                                                    control_col_ids, perc))
      outlist_noZ$sd.ctrls[metabolite_index]  <-   sd(robust_scaler(outlist_noZ[metabolite_index, control_col_ids],
                                                                    control_col_ids, perc))
    }
  }

  # Make and add columns with robust zscores
  cnames_robust <- gsub("_Zscore", "_RobustZscore", colnames_z)
  for (i in intensity_col_ids) {
    zscores_1col <- (as.numeric(as.vector(unlist(outlist_noZ[, i]))) - outlist_noZ$avg.ctrls) / outlist_noZ$sd.ctrls
    outlist_noZ <- cbind(outlist_noZ, zscores_1col)
  }
  colnames(outlist_noZ)[startcol:ncol(outlist_noZ)] <- cnames_robust

  # output metabolites filtered on relevance
  save(outlist, file = paste0("AdductSums_filtered_Zscores.RData"))
  write.table(outlist, file = paste0("AdductSums_filtered_Zscores.txt"), sep = "\t", row.names = FALSE)
  # output filtered metabolites with robust scaled Zscores
  save(outlist_noZ, file = paste0("AdductSums_filtered_robustZ.RData"))
  write.table(outlist_noZ, file = paste0("AdductSums_filtered_robustZ.txt"), sep = "\t", row.names = FALSE)

  # get the IDs of the patients and sort
  patient_ids <- unique(as.vector(unlist(lapply(strsplit(colnames(patient_columns), ".", fixed = TRUE), function(x) x[1]))))
  patient_ids <- patient_ids[order(nchar(patient_ids), patient_ids)]

  # for every row, make boxplot, insert into excel, and calculate Zscore for every patient
  temp_png <- NULL
  for (p in 1:nrow(outlist)) {
    # get HMDB ID
    hmdb_name <- rownames(outlist[p, ])
    # get intensities per metabolite for box plot for control samples
    intensities <- list(as.numeric(as.vector(unlist(control_columns[p, ]))))
    labels <- c("C", patient_ids)
    # get intensities per metabolite for box plot for patient samples
    for (i in 1:length(patient_ids)) {
      id <- patient_ids[i]
      # combine all intensities that start with the same string for patients
      # exception: if there is only one patient id, skip this step; nothing to combine
      if (ncol(patient_columns) > 1) {
        patient_int <- as.numeric(as.vector(unlist(outlist[p, names(patient_columns[1, ])
  						 [startsWith(names(patient_columns[1, ]), paste0(id, "."))]])))
      } else {
	patient_int <- as.numeric(unlist(as.vector(patient_columns)))
      }
      intensities[[i + 1]] <- patient_int
    }
    intensities <- setNames(intensities, labels)

    plot_width <- length(labels) * 16 + 90

    plot.new()
    if (export) {
      png(filename = paste0(plot_dir, "/", hmdb_name, "_box.png"),
          width = plot_width,
          height = 280)
    }
    # set margins
    par(oma = c(2, 0, 0, 0))
    boxplot(intensities,
            col = c("green", rep("red", length(intensities) - 1)),
            names.arg = labels,
            las = 2,
            main = hmdb_name)
    dev.off()

    file_png <- paste0(plot_dir, "/", hmdb_name, "_box.png")
    if (is.null(temp_png)) {
      temp_png <- loder::readPng(file_png)
      img_dim <- dim(temp_png)[c(1, 2)]
      cell_dim <- img_dim * imagesize_multiplier
      setColWidths(wb, filelist, cols = 1, widths = cell_dim[2] / 20)
    }

    openxlsx::insertImage(wb,
                          filelist,
                          file_png,
                          startRow = p + 1,
                          startCol = 1,
                          height = cell_dim[1],
                          width = cell_dim[2],
                          units = "px")

    if (p %% 100 == 0) {
      cat("at row: ", p, "\n")
    }
  }

  openxlsx::setRowHeights(wb, filelist, rows = c(1:nrow(outlist) + 1), heights = cell_dim[1] / 4)
  openxlsx::setColWidths(wb, filelist, cols = c(2:ncol(outlist)), widths = 20)
} else {
  openxlsx::setRowHeights(wb, filelist, rows = c(1:nrow(outlist)), heights = 18)
  openxlsx::setColWidths(wb, filelist, cols = c(1:ncol(outlist)), widths = 20)
}
print("Write Excel file")
# write Excel file
outlist$name <- stringi::stri_enc_toutf8(outlist$name)
openxlsx::writeData(wb, sheet = 1, outlist, startCol = 1)
xlsx_name <- paste0(outdir, "/", project, ".xlsx")
openxlsx::saveWorkbook(wb, xlsx_name, overwrite = TRUE)
rm(wb)

#### INTERNAL STANDARDS ####
is_list <- outlist[grep("Internal standard", outlist[, "relevance"], fixed = TRUE), ]
is_codes <- rownames(is_list)

# if all data from one samplename (for example P195.1) is filtered out in 3-averageTechReplicates 
# because of too little data (threshold parameter)i, the init.RData (repl_pattern) will contain more sample_names 
# than the peak data (IS), so this data needs to be removed first, before the retrieval of the summed adducts. 
# Write sample_names to a log file, to let user know that this sample_name contained no data.
sample_names_nodata <- setdiff(names(repl_pattern), names(is_list))
if (!is.null(sample_names_nodata)) {
  write.table(sample_names_nodata, file = paste(outdir, "sample_names_nodata.txt", sep = "/"), 
	      row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(sample_names_nodata, "\n")
  for (sample_name in sample_names_nodata) {
    repl_pattern[[sample_name]] <- NULL
  }
}

# Retrieve IS summed adducts
is_summed <- is_list[c(names(repl_pattern), "HMDB_code")]
is_summed$HMDB.name <- is_list$name
is_summed <- reshape2::melt(is_summed, id.vars = c("HMDB_code", "HMDB.name"))
colnames(is_summed) <- c("HMDB.code", "HMDB.name", "Sample", "Intensity")
is_summed$Matrix <- dims_matrix
is_summed$Rundate <- rundate
is_summed$Project <- project
is_summed$Intensity <- as.numeric(as.character(is_summed$Intensity))

# Retrieve IS positive mode
is_pos <- as.data.frame(subset(outlist_pos_adducts_hmdb, rownames(outlist_pos_adducts_hmdb) %in% is_codes))
is_pos <- is_pos[c(names(repl_pattern))]
is_pos$HMDB_name <- is_list[match(row.names(is_pos), is_list$HMDB_code, nomatch = NA), "name"]
is_pos$HMDB.code <- row.names(is_pos)
is_pos <- reshape2::melt(is_pos, id.vars = c("HMDB.code", "HMDB_name"))
colnames(is_pos) <- c("HMDB.code", "HMDB.name", "Sample", "Intensity")
is_pos$Matrix <- dims_matrix
is_pos$Rundate <- rundate
is_pos$Project <- project
is_pos$Intensity <- as.numeric(as.character(is_pos$Intensity))

# Retrieve IS negative mode
is_neg <- as.data.frame(subset(outlist_neg_adducts_hmdb, rownames(outlist_neg_adducts_hmdb) %in% is_codes))
is_neg <- is_neg[c(names(repl_pattern))]
is_neg$HMDB_name <- is_list[match(row.names(is_neg), is_list$HMDB_code, nomatch = NA), "name"]
is_neg$HMDB.code <- row.names(is_neg)
is_neg <- reshape2::melt(is_neg, id.vars = c("HMDB.code", "HMDB_name"))
colnames(is_neg) <- c("HMDB.code", "HMDB.name", "Sample", "Intensity")
is_neg$Matrix <- dims_matrix
is_neg$Rundate <- rundate
is_neg$Project <- project
is_neg$Intensity <- as.numeric(as.character(is_neg$Intensity))
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
# theme for all IS bar plots
theme_is_bar <- function(my_plot) {
  my_plot +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::theme(legend.position = "none", 
	           axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), 
               	   axis.text.y = element_text(size = 6)
    )
}

# ggplot functions
is_neg_bar_plot <- ggplot(is_neg, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Neg)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free_y")

is_pos_bar_plot <- ggplot(is_pos, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Pos)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free_y")

is_sum_bar_plot <- ggplot(is_summed, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Summed)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free_y")

# add theme to ggplot functions
is_neg_bar_plot <- theme_is_bar(is_neg_bar_plot)
is_pos_bar_plot <- theme_is_bar(is_pos_bar_plot)
is_sum_bar_plot <- theme_is_bar(is_sum_bar_plot)

# save plots to disk
plot_width <- 9 + 0.35 * sample_count
ggsave(paste0(outdir, "/plots/IS_bar_all_neg.png"), 
       plot = is_neg_bar_plot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_bar_all_pos.png"), 
       plot = is_pos_bar_plot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_bar_all_sum.png"), 
       plot = is_sum_bar_plot, height = plot_width / 2.5, width = plot_width, units = "in")

## Line plots with all IS
# function for ggplot theme
# add smaller legend in the "all IS line plots", otherwise out-of-range when more than 13 IS lines
theme_is_line <- function(my_plot) {
  my_plot +
    guides(shape = guide_legend(override.aes = list(size = 0.5)), 
	   color = guide_legend(override.aes = list(size = 0.5))
    ) +
    theme(legend.title = element_text(size = 8), 
	  legend.text = element_text(size = 6), 
	  legend.key.size = unit(0.7, "line"), 
	  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)
    )
}

# ggplot functions
is_neg_line_plot <- ggplot(is_neg, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Neg)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

is_pos_line_plot <- ggplot(is_pos, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Pos)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

is_sum_line_plot <- ggplot(is_summed, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Sum)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

# add theme to ggplot functions
is_sum_line_plot <- theme_is_line(is_sum_line_plot)
is_neg_line_plot <- theme_is_line(is_neg_line_plot)
is_pos_line_plot <- theme_is_line(is_pos_line_plot)

# save plots to disk
plot_width <- 8 + 0.2 * sample_count
ggsave(paste0(outdir, "/plots/IS_line_all_neg.png"), 
       plot = is_neg_line_plot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_line_all_pos.png"), 
       plot = is_pos_line_plot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_line_all_sum.png"), 
       plot = is_sum_line_plot, height = plot_width / 2.5, width = plot_width, units = "in")

## bar plots with a selection of IS
is_neg_selection <- c("2H2-Ornithine (IS)", "2H3-Glutamate (IS)", "2H2-Citrulline (IS)", "2H4_13C5-Arginine (IS)", 
		      "13C6-Tyrosine (IS)")
is_pos_selection <- c("2H4-Alanine (IS)", "13C6-Phenylalanine (IS)", "2H4_13C5-Arginine (IS)", "2H3-Propionylcarnitine (IS)", 
		      "2H9-Isovalerylcarnitine (IS)")
is_sum_selection <- c("2H8-Valine (IS)", "2H3-Leucine (IS)", "2H3-Glutamate (IS)", "2H4_13C5-Arginine (IS)", 
		      "13C6-Tyrosine (IS)")

# add minimal intensity lines based on matrix (DBS or Plasma) and machine mode (neg, pos, sum)
if (dims_matrix == "DBS") {
  hline_data_neg <-
    data.frame(
      z = c(15000, 200000, 130000, 18000, 50000),
      HMDB.name = is_neg_selection
    )
  hline_data_pos <-
    data.frame(
      z = c(150000, 3300000, 1750000, 150000, 270000),
      HMDB.name = is_pos_selection
    )
  hline_data_sum <-
    data.frame(
      z = c(1300000, 2500000, 500000, 1800000, 1400000),
      HMDB.name = is_sum_selection
    )
} else if (dims_matrix == "Plasma") {
  hline_data_neg <-
    data.frame(
      z = c(6500, 100000, 75000, 7500, 25000),
      HMDB.name = is_neg_selection
    )
  hline_data_pos <-
    data.frame(
      z = c(85000, 1000000, 425000, 70000, 180000),
      HMDB.name = is_pos_selection
    )
  hline_data_sum <-
    data.frame(
      z = c(700000, 1250000, 150000, 425000, 300000),
      HMDB.name = is_sum_selection
    )
}

# ggplot functions
is_neg_selection_barplot <- ggplot(subset(is_neg, HMDB.name %in% is_neg_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Neg)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free", ncol = 2) +
  if (exists("hline_data_neg")) {
    geom_hline(aes(yintercept = z), subset(hline_data_neg, HMDB.name %in% is_neg$HMDB.name))
  }

is_pos_selection_barplot <- ggplot(subset(is_pos, HMDB.name %in% is_pos_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Pos)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free", ncol = 2) +
  if (exists("hline_data_pos")) {
    geom_hline(aes(yintercept = z), subset(hline_data_pos, HMDB.name %in% is_pos$HMDB.name))
  }

is_sum_selection_barplot <- ggplot(subset(is_summed, HMDB.name %in% is_sum_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Sum)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free", ncol = 2) +
  if (exists("hline_data_sum")) {
    geom_hline(aes(yintercept = z), subset(hline_data_sum, HMDB.name %in% is_summed$HMDB.name))
  }

# add theme to ggplot functions
is_neg_selection_barplot <- theme_is_bar(is_neg_selection_barplot)
is_pos_selection_barplot <- theme_is_bar(is_pos_selection_barplot)
is_sum_selection_barplot <- theme_is_bar(is_sum_selection_barplot)

# save plots to disk
plot_width <- 9 + 0.35 * sample_count
ggsave(paste0(outdir, "/plots/IS_bar_select_neg.png"), 
       plot = is_neg_selection_barplot, height = plot_width / 2.0, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_bar_select_pos.png"), 
       plot = is_pos_selection_barplot, height = plot_width / 2.0, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_bar_select_sum.png"), 
       plot = is_sum_selection_barplot, height = plot_width / 2.0, width = plot_width, units = "in")

## line plots with a selection of IS
# ggplot functions
is_neg_selection_lineplot <- ggplot(subset(is_neg, HMDB.name %in% is_neg_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Neg)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

is_pos_selection_lineplot <- ggplot(subset(is_pos, HMDB.name %in% is_pos_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Pos)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

is_sum_selection_lineplot <- ggplot(subset(is_summed, HMDB.name %in% is_sum_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Sum)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

# add theme to ggplot functions
is_neg_selection_lineplot <- theme_is_line(is_neg_selection_lineplot)
is_pos_selection_lineplot <- theme_is_line(is_pos_selection_lineplot)
is_sum_selection_lineplot <- theme_is_line(is_sum_selection_lineplot)

# save plots to disk
plot_width <- 8 + 0.2 * sample_count
ggsave(paste0(outdir, "/plots/IS_line_select_neg.png"), 
       plot = is_neg_selection_lineplot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_line_select_pos.png"), 
       plot = is_pos_selection_lineplot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_line_select_sum.png"), 
       plot = is_sum_selection_lineplot, height = plot_width / 2.5, width = plot_width, units = "in")


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
  if (any(grep("^(P1002\\.)[[:digit:]]+_", positive_control_list)) &&
    any(grep("^(P1003\\.)[[:digit:]]+_", positive_control_list)) &&
    any(grep("^(P1005\\.)[[:digit:]]+_", positive_control_list))) {
    cat("All three positive controls are present")
  } else {
    pos_contr_warning <- paste0(c("positive controls list is not complete. Only ", 
				  positive_control_list, " is/are present"), collapse = " ")
  }
  # you need all positive control samples, thus starting the script only if all are available
  if (length(pos_contr_warning) == 0) {
    # make positive control excel with specific HMDB_codes in combination with specific control samples
    pa_sample_name <- positive_control_list[grepl("P1002", positive_control_list)]
    pku_sample_name <- positive_control_list[grepl("P1003", positive_control_list)]
    lpi_sample_name <- positive_control_list[grepl("P1005", positive_control_list)]

    # pa_codes <- c("HMDB0000824", "HMDB0000783", "HMDB0000123")
    # pku_codes <- c("HMDB0000159")
    # lpi_codes <- c("HMDB0000904", "HMDB0000641", "HMDB0000182")

    pa_codes <- c("HMDB0000824", "HMDB0000725", "HMDB0000123")
    pku_codes <- c("HMDB0000159")
    lpi_codes <- c("HMDB0000904", "HMDB0000641", "HMDB0000182")

    pa_names <- c("Propionylcarnitine", "Propionylglycine", "Glycine")
    pku_names <- c("L-Phenylalanine")
    lpi_names <- c("Citrulline", "L-Glutamine", "L-Lysine")

    pa_data <- outlist[pa_codes, c("HMDB_code", "name", pa_sample_name)]
    pa_data <- reshape2::melt(pa_data, id.vars = c("HMDB_code", "name"))
    colnames(pa_data) <- c("HMDB.code", "HMDB.name", "Sample", "Zscore")
    # change HMDB names because propionylglycine doesn't have its own line, rowname is HMDB0000725 (4-hydroxyproline)
    pa_data$HMDB.name <- pa_codes
    # change HMDB codes so the propionylglycine is the correct HMDB ID
    pa_data$HMDB.code <- c("HMDB0000824", "HMDB0000783", "HMDB0000123")

    pku_data <- outlist[pku_codes, c("HMDB_code", "name", pku_sample_name)]
    pku_data <- reshape2::melt(pku_data, id.vars = c("HMDB_code", "name"))
    colnames(pku_data) <- c("HMDB.code", "HMDB.name", "Sample", "Zscore")

    lpi_data <- outlist[lpi_codes, c("HMDB_code", "name", lpi_sample_name)]
    lpi_data <- reshape2::melt(lpi_data, id.vars = c("HMDB_code", "name"))
    colnames(lpi_data) <- c("HMDB.code", "HMDB.name", "Sample", "Zscore")

    positive_control <- rbind(pa_data, pku_data, lpi_data)
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
	       sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)
  } else {
    write.table(pos_contr_warning, file = paste(outdir, "positive_controls_warning.txt", sep = "/"), 
		row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

### MISSING M/Z CHECK
# check the outlist_identified_(negative/positive).RData files for missing m/z values and mention in the results mail
# Load the outlist_identified files + remove the loaded files
load(paste0(outdir, "/outlist_identified_negative.RData"))
outlist_ident_neg <- outlist_ident
load(paste0(outdir, "/outlist_identified_positive.RData"))
outlist_ident_pos <- outlist_ident
rm(outlist_ident)
# check for missing m/z in negative and positive mode
scanmode <- c("Negative", "Positive")
index <- 1
results_ident <- c() 
outlist_ident_list <- list(outlist_ident_neg, outlist_ident_pos)
for (outlist_ident in outlist_ident_list) {
  current_mode <- scanmode[index]
  # retrieve all unique m/z values in whole numbers and check if all are available
  mz_values <- as.numeric(unique(format(outlist_ident$mzmed.pgrp, digits = 0)))
  # m/z range for a standard run = 70-600
  mz_range <- seq(70, 599, by = 1)
  mz_missing <- c()
  for (mz in mz_range) {
    if (!mz %in% mz_values) {
      mz_missing <- c(mz_missing, mz)
    }
  }
  y <- mz_missing
  # check if m/z are missing and make an .txt file with information
  group_ident <- cumsum(c(1, abs(y[-length(y)] - y[-1]) > 1))
  if (length(group_ident) > 1) {
    results_ident <- c(results_ident, paste0("Missing m/z values ", current_mode, " mode"))
    results_ident <- c(results_ident, by(y, group_ident, identity))
  } else {
    results_ident <- c(results_ident, paste0(current_mode, " mode did not have missing mz values"))
  }
  # change to other scanmode
  index <- index + 1
}
lapply(results_ident, write, file = paste0(outdir, "/missing_mz_warning.txt"), append = TRUE, ncolumns = 1000)
