#!/usr/bin/Rscript
## adapted from 13-excelExport.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

init_filepath <- cmd_args[1]
project       <- cmd_args[2]
dims_matrix   <- cmd_args[3]
hmdb          <- cmd_args[4]
z_score       <- as.numeric(cmd_args[5])

# load required packages
library("ggplot2")
library("reshape2")
library("openxlsx")
library("loder")
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))

round_df <- function(df, digits) {
  #' function for rounding numbers to x digits for numeric values
  #'
  #' @param df dataframe
  #' @param digits integer number of digits to round off to
  #'
  #' @return df
  numeric_columns <- sapply(df, mode) == "numeric"
  df[numeric_columns] <- round(df[numeric_columns], digits)
  return(df)
}

robust_scaler <- function(control_intensities, control_col_ids, perc = 5) {
  #' calculate robust scaler: Z-score based on controls without outliers
  #'
  #' @param control_intensities matrix with intensities for control samples
  #' @param control_col_ids vector with column names for control samples
  #' @param perc float percentage of outliers which will be removed from controls
  #'
  #' @return trimmed_control_intensities
  nr_toremove <- ceiling(length(control_col_ids) * perc / 100)
  sorted_control_intensities <- sort(as.numeric(control_intensities))
  trimmed_control_intensities <- sorted_control_intensities[(nr_toremove + 1) : (length(sorted_control_intensities) - nr_toremove)]
  return(trimmed_control_intensities)
}

# Initialise and load data
plot   <- TRUE
export <- TRUE
control_label <- "C"
case_label    <- "P"
imagesize_multiplier <- 2
# setting outdir to export files to the project directory
outdir <- "./"
# percentage of outliers to remove from calculation of robust scaler
perc <- 5

# load information on samples
load(init_filepath)
# load the HMDB file with info on biological relevance of metabolites
load(hmdb)

# get current date
rundate <- Sys.Date()

# create a directory for plots in project directory
plotdir <- paste0(outdir, "/plots/adducts")
dir.create(paste0(outdir, "/plots"), showWarnings = FALSE)
dir.create(plotdir, showWarnings = FALSE)

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
tmp_pos <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], 1:(dim(outlist_pos_adducts_hmdb)[2] - 1)]
tmp_hmdb_name_pos <- outlist_pos_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], dim(outlist_pos_adducts_hmdb)[2]]
tmp_pos_left <- outlist_pos_adducts_hmdb[-index_pos, ]
# same for negative mode
tmp_neg <- outlist_neg_adducts_hmdb[rownames(outlist_pos_adducts_hmdb)[index_pos], 1:(dim(outlist_neg_adducts_hmdb)[2] - 1)]
tmp_neg_left <- outlist_neg_adducts_hmdb[-index_neg, ]

# Combine positive and negative numbers and paste back HMDB column
tmp <- apply(tmp_pos, 2, as.numeric) + apply(tmp_neg, 2, as.numeric)
rownames(tmp) <- rownames(tmp_pos)
tmp <- cbind(tmp, "HMDB_name" = tmp_hmdb_name_pos)
outlist <- rbind(tmp, tmp_pos_left, tmp_neg_left)

# Filter for biological relevance
peaks_in_list <- which(rownames(outlist) %in% rownames(rlvnc))
outlist <- cbind(outlist[peaks_in_list, ], as.data.frame(rlvnc[rownames(outlist)[peaks_in_list], ]))
# filter out all irrelevant HMDBs
outlist <- outlist %>%
  tibble::rownames_to_column("rowname") %>%
  filter(!grepl("Exogenous|Drug|exogenous", relevance)) %>%
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
  control_columns <- outlist[, control_col_ids]

  # Get columns with patient intensities
  patient_col_ids <- grep(case_label, colnames(outlist), fixed = TRUE)
  patient_columns <- outlist[, patient_col_ids]
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

  # output metabolites filtered on relevance into tab-separated file
  write.table(outlist, file = paste0("AdductSums_filtered_Zscores.txt"), sep = "\t", row.names = FALSE)
  # output filtered metabolites with robust scaled Zscores
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
      patient_int <- as.numeric(as.vector(unlist(outlist[p, names(patient_columns[1, ])[startsWith(names(patient_columns[1, ]), paste0(id, "."))]])))
      intensities[[i + 1]] <- patient_int
    }
    intensities <- setNames(intensities, labels)

    plot_width <- length(labels) * 16 + 90

    plot.new()
    if (export) {
      png(filename = paste0(plotdir, "/", hmdb_name, "_box.png"),
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

    file_png <- paste0(plotdir, "/", hmdb_name, "_box.png")
    if (is.null(temp_png)) {
      temp_png <- readPng(file_png)
      img_dim <- dim(temp_png)[c(1, 2)]
      cell_dim <- img_dim * imagesize_multiplier
      setColWidths(wb, filelist, cols = 1, widths = cell_dim[2] / 20)
    }

    insertImage(wb,
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

  setRowHeights(wb, filelist, rows = c(1:nrow(outlist) + 1), heights = cell_dim[1] / 4)
  setColWidths(wb, filelist, cols = c(2:ncol(outlist)), widths = 20)
} else {
  setRowHeights(wb, filelist, rows = c(1:nrow(outlist)), heights = 18)
  setColWidths(wb, filelist, cols = c(1:ncol(outlist)), widths = 20)
}
# write Excel file
writeData(wb, sheet = 1, outlist, startCol = 1)
xlsx_name <- paste0(outdir, "/", project, ".xlsx")
saveWorkbook(wb, xlsx_name, overwrite = TRUE)
cat(xlsx_name)
rm(wb)

#### INTERNE STANDAARDEN ####
IS <- outlist[grep("Internal standard", outlist[, "relevance"], fixed = TRUE), ]
IS_codes <- rownames(IS)
cat(IS_codes, "\n")

# if all data from one samplename (for example P195.1) is filtered out in 3-averageTechReplicates because of too little data (threshold parameter) the init.RData (repl_pattern) will contain more sample_names then the peak data (IS),
# thus this data needs to be removed first, before the retrieval of the summed adducts. Write sample_names to a log file, to let user know that this sample_name contained no data.
sample_names_nodata <- setdiff(names(repl_pattern), names(IS))
if (!is.null(sample_names_nodata)) {
  write.table(sample_names_nodata, file = paste(outdir, "sample_names_nodata.txt", sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  cat(sample_names_nodata, "\n")
  for (sample_name in sample_names_nodata) {
    repl_pattern[[sample_name]] <- NULL
  }
}

# Retrieve IS summed adducts
IS_summed <- IS[c(names(repl_pattern), "HMDB_code")]
IS_summed$HMDB.name <- IS$name
IS_summed <- reshape2::melt(IS_summed, id.vars = c("HMDB_code", "HMDB.name"))
colnames(IS_summed) <- c("HMDB.code", "HMDB.name", "Sample", "Intensity")
IS_summed$Intensity <- as.numeric(IS_summed$Intensity)
IS_summed$Matrix <- dims_matrix
IS_summed$Rundate <- rundate
IS_summed$Project <- project
IS_summed$Intensity <- as.numeric(as.character(IS_summed$Intensity))

# Retrieve IS positive mode
IS_pos <- as.data.frame(subset(outlist_pos_adducts_hmdb, rownames(outlist_pos_adducts_hmdb) %in% IS_codes))
IS_pos$HMDB_name <- IS[match(row.names(IS_pos), IS$HMDB_code, nomatch = NA), "name"]
IS_pos$HMDB.code <- row.names(IS_pos)
IS_pos <- reshape2::melt(IS_pos, id.vars = c("HMDB.code", "HMDB_name"))
colnames(IS_pos) <- c("HMDB.code", "HMDB.name", "Sample", "Intensity")
IS_pos$Matrix <- dims_matrix
IS_pos$Rundate <- rundate
IS_pos$Project <- project
IS_pos$Intensity <- as.numeric(as.character(IS_pos$Intensity))

# Retrieve IS negative mode
IS_neg <- as.data.frame(subset(outlist_neg_adducts_hmdb, rownames(outlist_neg_adducts_hmdb) %in% IS_codes))
IS_neg$HMDB_name <- IS[match(row.names(IS_neg), IS$HMDB_code, nomatch = NA), "name"]
IS_neg$HMDB.code <- row.names(IS_neg)
IS_neg <- reshape2::melt(IS_neg, id.vars = c("HMDB.code", "HMDB_name"))
colnames(IS_neg) <- c("HMDB.code", "HMDB.name", "Sample", "Intensity")
IS_neg$Matrix <- dims_matrix
IS_neg$Rundate <- rundate
IS_neg$Project <- project
IS_neg$Intensity <- as.numeric(as.character(IS_neg$Intensity))

# Save results
save(IS_pos, IS_neg, IS_summed, file = paste0(outdir, "/", project, "_IS_results.RData"))

# number of samples, for plotting length and width
sample_count <- length(repl_pattern)

# change the order of the x-axis summed plots to a natural sorted one
Sample_naturalorder <- unique(as.character(IS_summed$Sample))
Sample_naturalorder <- str_sort(Sample_naturalorder, numeric = TRUE)
IS_summed$Sample_level <- factor(IS_summed$Sample, levels = c(Sample_naturalorder))
IS_pos$Sample_level <- factor(IS_pos$Sample, levels = c(Sample_naturalorder))
IS_neg$Sample_level <- factor(IS_neg$Sample, levels = c(Sample_naturalorder))

## bar plots with all IS

# function for ggplot theme
# theme for all IS bar plots
theme_IS_bar <- function(myPlot) {
  myPlot +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y = element_text(size = 6)
    )
}

# ggplot functions
IS_neg_bar_plot <- ggplot(IS_neg, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Neg)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free_y")

IS_pos_bar_plot <- ggplot(IS_pos, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Pos)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free_y")

IS_sum_bar_plot <- ggplot(IS_summed, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Summed)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free_y")

# add theme to ggplot functions
IS_neg_bar_plot <- theme_IS_bar(IS_neg_bar_plot)
IS_pos_bar_plot <- theme_IS_bar(IS_pos_bar_plot)
IS_sum_bar_plot <- theme_IS_bar(IS_sum_bar_plot)

# save plots to disk
plot_width <- 9 + 0.35 * sample_count
ggsave(paste0(outdir, "/plots/IS_bar_all_neg.png"), plot = IS_neg_bar_plot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_bar_all_pos.png"), plot = IS_pos_bar_plot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_bar_all_sum.png"), plot = IS_sum_bar_plot, height = plot_width / 2.5, width = plot_width, units = "in")

## Line plots with all IS

# function for ggplot theme
# add smaller legend in the "all IS line plots", otherwise out-of-range when more than 13 IS lines
theme_IS_line <- function(myPlot) {
  myPlot +
    guides(
      shape = guide_legend(override.aes = list(size = 0.5)),
      color = guide_legend(override.aes = list(size = 0.5))
    ) +
    theme(
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.7, "line"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)
    )
}

# ggplot functions
IS_neg_line_plot <- ggplot(IS_neg, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Neg)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

IS_pos_line_plot <- ggplot(IS_pos, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Pos)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

IS_sum_line_plot <- ggplot(IS_summed, aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Sum)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

# add theme to ggplot functions
IS_sum_line_plot <- theme_IS_line(IS_sum_line_plot)
IS_neg_line_plot <- theme_IS_line(IS_neg_line_plot)
IS_pos_line_plot <- theme_IS_line(IS_pos_line_plot)

# save plots to disk
plot_width <- 8 + 0.2 * sample_count
ggsave(paste0(outdir, "/plots/IS_line_all_neg.png"), plot = IS_neg_line_plot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_line_all_pos.png"), plot = IS_pos_line_plot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_line_all_sum.png"), plot = IS_sum_line_plot, height = plot_width / 2.5, width = plot_width, units = "in")

## bar plots with a selection of IS
IS_neg_selection <- c("2H2-Ornithine (IS)", "2H3-Glutamate (IS)", "2H2-Citrulline (IS)", "2H4_13C5-Arginine (IS)", "13C6-Tyrosine (IS)")
IS_pos_selection <- c("2H4-Alanine (IS)", "13C6-Phenylalanine (IS)", "2H4_13C5-Arginine (IS)", "2H3-Propionylcarnitine (IS)", "2H9-Isovalerylcarnitine (IS)")
IS_sum_selection <- c("2H8-Valine (IS)", "2H3-Leucine (IS)", "2H3-Glutamate (IS)", "2H4_13C5-Arginine (IS)", "13C6-Tyrosine (IS)")

# add minimal intensity lines based on matrix (DBS or Plasma) and machine mode (neg, pos, sum)
if (dims_matrix == "DBS") {
  hline.data.neg <-
    data.frame(
      z = c(15000, 200000, 130000, 18000, 50000),
      HMDB.name = IS_neg_selection
    )
  hline.data.pos <-
    data.frame(
      z = c(150000, 3300000, 1750000, 150000, 270000),
      HMDB.name = IS_pos_selection
    )
  hline.data.sum <-
    data.frame(
      z = c(1300000, 2500000, 500000, 1800000, 1400000),
      HMDB.name = IS_sum_selection
    )
} else if (dims_matrix == "Plasma") {
  hline.data.neg <-
    data.frame(
      z = c(6500, 100000, 75000, 7500, 25000),
      HMDB.name = IS_neg_selection
    )
  hline.data.pos <-
    data.frame(
      z = c(85000, 1000000, 425000, 70000, 180000),
      HMDB.name = IS_pos_selection
    )
  hline.data.sum <-
    data.frame(
      z = c(700000, 1250000, 150000, 425000, 300000),
      HMDB.name = IS_sum_selection
    )
}

# function for ggplot theme
# see bar plots with all IS

# ggplot functions
IS_neg_selection_barplot <- ggplot(subset(IS_neg, HMDB.name %in% IS_neg_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Neg)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free", ncol = 2) +
  if (exists("hline.data.neg")) {
    geom_hline(aes(yintercept = z), subset(hline.data.neg, HMDB.name %in% IS_neg$HMDB.name))
  } # subset, if some IS have no data, no empty plots will be generated with a line)

IS_pos_selection_barplot <- ggplot(subset(IS_pos, HMDB.name %in% IS_pos_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Pos)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free", ncol = 2) +
  if (exists("hline.data.pos")) {
    geom_hline(aes(yintercept = z), subset(hline.data.pos, HMDB.name %in% IS_pos$HMDB.name))
  }

IS_sum_selection_barplot <- ggplot(subset(IS_summed, HMDB.name %in% IS_sum_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Sum)") +
  geom_bar(aes(fill = HMDB.name), stat = "identity") +
  labs(x = "", y = "Intensity") +
  facet_wrap(~HMDB.name, scales = "free", ncol = 2) +
  if (exists("hline.data.sum")) {
    geom_hline(aes(yintercept = z), subset(hline.data.sum, HMDB.name %in% IS_summed$HMDB.name))
  }

# add theme to ggplot functions
IS_neg_selection_barplot <- theme_IS_bar(IS_neg_selection_barplot)
IS_pos_selection_barplot <- theme_IS_bar(IS_pos_selection_barplot)
IS_sum_selection_barplot <- theme_IS_bar(IS_sum_selection_barplot)

# save plots to disk
plot_width <- 9 + 0.35 * sample_count
ggsave(paste0(outdir, "/plots/IS_bar_select_neg.png"), plot = IS_neg_selection_barplot, height = plot_width / 2.0, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_bar_select_pos.png"), plot = IS_pos_selection_barplot, height = plot_width / 2.0, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_bar_select_sum.png"), plot = IS_sum_selection_barplot, height = plot_width / 2.0, width = plot_width, units = "in")

## line plots with a selection of IS

# function for ggplot theme
# see line plots with all IS

# ggplot functions
IS_neg_selection_lineplot <- ggplot(subset(IS_neg, HMDB.name %in% IS_neg_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Neg)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

IS_pos_selection_lineplot <- ggplot(subset(IS_pos, HMDB.name %in% IS_pos_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Pos)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

IS_sum_selection_lineplot <- ggplot(subset(IS_summed, HMDB.name %in% IS_sum_selection), aes(Sample_level, Intensity)) +
  ggtitle("Interne Standaard (Sum)") +
  geom_point(aes(col = HMDB.name)) +
  geom_line(aes(col = HMDB.name, group = HMDB.name)) +
  labs(x = "", y = "Intensity")

# add theme to ggplot functions
IS_neg_selection_lineplot <- theme_IS_line(IS_neg_selection_lineplot)
IS_pos_selection_lineplot <- theme_IS_line(IS_pos_selection_lineplot)
IS_sum_selection_lineplot <- theme_IS_line(IS_sum_selection_lineplot)

# save plots to disk
plot_width <- 8 + 0.2 * sample_count
ggsave(paste0(outdir, "/plots/IS_line_select_neg.png"), plot = IS_neg_selection_lineplot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_line_select_pos.png"), plot = IS_pos_selection_lineplot, height = plot_width / 2.5, width = plot_width, units = "in")
ggsave(paste0(outdir, "/plots/IS_line_select_sum.png"), plot = IS_sum_selection_lineplot, height = plot_width / 2.5, width = plot_width, units = "in")


### POSITIVE CONTROLS CHECK
# these positive controls need to be in the samplesheet, in order to make the Pos_Contr.RData file
# Positive control samples all have the format P1002.x, P1003.x and P1005.x (where x is a number)

column_list <- colnames(outlist)
patterns <- c("^(P1002\\.)[[:digit:]]+_", "^(P1003\\.)[[:digit:]]+_", "^(P1005\\.)[[:digit:]]+_")
positive_controls_index <- grepl(pattern = paste(patterns, collapse = "|"), column_list)
positivecontrol_list <- column_list[positive_controls_index]

if (z_score == 1) {
  # find if one or more positive control samples are missing
  pos_contr_warning <- c()
  # any() grep because you get a vector of FALSE's and TRUE's. only one grep match is needed for each positive control
  if (any(grep("^(P1002\\.)[[:digit:]]+_", positivecontrol_list)) &&
    any(grep("^(P1003\\.)[[:digit:]]+_", positivecontrol_list)) &&
    any(grep("^(P1005\\.)[[:digit:]]+_", positivecontrol_list))) {
    cat("All three positive controls are present")
  } else {
    pos_contr_warning <- paste0(c("positive controls list is not complete. Only ", positivecontrol_list, " is/are present"), collapse = " ")
  }
  # you need all positive control samples, thus starting the script only if all are available
  if (length(pos_contr_warning) == 0) {
    ### POSITIVE CONTROLS
    # make positive control excel with specific HMDB_codes in combination with specific control samples
    PA_sample_name <- positivecontrol_list[grepl("P1002", positivecontrol_list)] # P1001.x_Zscore
    PKU_sample_name <- positivecontrol_list[grepl("P1003", positivecontrol_list)] # P1003.x_Zscore
    LPI_sample_name <- positivecontrol_list[grepl("P1005", positivecontrol_list)] # P1005.x_Zscore

    PA_codes <- c("HMDB00824", "HMDB00783", "HMDB00123")
    PKU_codes <- c("HMDB00159")
    LPI_codes <- c("HMDB00904", "HMDB00641", "HMDB00182")

    PA_data <- outlist[PA_codes, c("HMDB_code", "name", PA_sample_name)]
    PA_data <- reshape2::melt(PA_data, id.vars = c("HMDB_code", "name"))
    colnames(PA_data) <- c("HMDB.code", "HMDB.name", "Sample", "Zscore")

    PKU_data <- outlist[PKU_codes, c("HMDB_code", "name", PKU_sample_name)]
    PKU_data <- reshape2::melt(PKU_data, id.vars = c("HMDB_code", "name"))
    colnames(PKU_data) <- c("HMDB.code", "HMDB.name", "Sample", "Zscore")

    LPI_data <- outlist[LPI_codes, c("HMDB_code", "name", LPI_sample_name)]
    LPI_data <- reshape2::melt(LPI_data, id.vars = c("HMDB_code", "name"))
    colnames(LPI_data) <- c("HMDB.code", "HMDB.name", "Sample", "Zscore")

    Pos_Contr <- rbind(PA_data, PKU_data, LPI_data)
    Pos_Contr$Zscore <- as.numeric(Pos_Contr$Zscore)
    # extra information added to excel for future reference. made in beginning of this script
    Pos_Contr$Matrix <- dims_matrix
    Pos_Contr$Rundate <- rundate
    Pos_Contr$Project <- project

    # Save results
    save(Pos_Contr, file = paste0(outdir, "/", project, "_Pos_Contr.RData"))
    Pos_Contr$Zscore <- round_df(Pos_Contr$Zscore, 2) # asked by Lab to round the number to 2 digits
    write.xlsx(Pos_Contr, file = paste0(outdir, "/", project, "_Pos_Contr.xlsx"), sheetName = "Sheet1", col.names = TRUE, row.names = TRUE, append = FALSE)
  } else {
    write.table(pos_contr_warning, file = paste(outdir, "positive_controls_warning.txt", sep = "/"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

### MISSING M/Z CHECK
# check the outlist_identified_(negative/positive).RData files for missing m/z values and mention in the results mail
print("Nu in missing m/z check")
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
