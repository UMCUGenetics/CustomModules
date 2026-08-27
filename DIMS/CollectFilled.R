# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

preprocessing_scripts_dir <- cmd_args[1]
z_score <- as.numeric(cmd_args[2])

source(paste0(preprocessing_scripts_dir, "collect_filled_functions.R"))

# Initialize
options(digits = 16)

# for each scan mode, collect all filled peak group lists
scanmodes <- c("positive", "negative")
for (scanmode in scanmodes) {
  # get list of files
  filled_files <- list.files("./", full.names = TRUE, pattern = paste0(scanmode, "_identified_filled"))
  # load files and combine into one object
  peakgroup_list_total <- NULL
  for (file_nr in seq_along(filled_files)) {
    peakgrouplist_filled <- get(load(filled_files[file_nr]))
    peakgroup_list_total <- rbind(peakgroup_list_total, peakgrouplist_filled)
  }
  # remove duplicates; peak groups with exactly the same m/z
  peakgroup_list_total <- merge_duplicate_rows(peakgroup_list_total)
  # sort on mass
  peakgroup_list_total <- peakgroup_list_total[order(peakgroup_list_total[, "mzmed.pgrp"]), ]
  # load replication pattern
  pattern_file <- paste0(scanmode, "_repl_pattern.RData")
  repl_pattern <- get(load(pattern_file))
  # calculate Z-scores
  if (z_score == 1) {
    peakgroup_list_stats <- calculate_zscores_peakgrouplist(peakgroup_list_total)
  } else {
    peakgroup_list_stats <- peakgroup_list_total
  }
  # calculate ppm deviation
  peakgroup_list_withppm <- calculate_ppm_deviation(peakgroup_list_stats)
  #  put columns in correct order
  peakgroup_list_ident <- order_columns_peakgrouplist(peakgroup_list_withppm)

  # generate output in Excel-readable format:
  remove_columns <- c("mzmin.pgrp", "mzmax.pgrp")
  peakgroup_list_ident <- peakgroup_list_ident[, -which(colnames(peakgroup_list_ident) %in% remove_columns)]
  write.table(peakgroup_list_ident, file = paste0("peakgroup_list_identified_", scanmode, ".txt"), sep = "\t", row.names = FALSE)
  # export output in RData format
  save(peakgroup_list_ident, file = paste0("peakgroup_list_identified_", scanmode, ".RData"))
}
