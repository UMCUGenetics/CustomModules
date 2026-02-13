# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

preprocessing_scripts_dir <- cmd_args[1]
ppm <- as.numeric(cmd_args[2])
z_score <- as.numeric(cmd_args[3])

source(paste0(preprocessing_scripts_dir, "collect_filled_functions.R"))

# for each scan mode, collect all filled peak group lists
scanmodes <- c("positive", "negative")
for (scanmode in scanmodes) {
  # get list of files
  filled_files <- list.files("./", full.names = TRUE, pattern = paste0(scanmode, "_identified_filled"))
  # load files and combine into one object
  outlist_total <- NULL
  for (file_nr in seq_along(filled_files)) {
    peakgrouplist_filled <- get(load(filled_files[file_nr]))
    outlist_total <- rbind(outlist_total, peakgrouplist_filled)
  }
  # remove duplicates; peak groups with exactly the same m/z
  outlist_total <- merge_duplicate_rows(outlist_total)
  # sort on mass
  outlist_total <- outlist_total[order(outlist_total[, "mzmed.pgrp"]), ]
  # load replication pattern
  pattern_file <- paste0(scanmode, "_repl_pattern.RData")
  repl_pattern <- get(load(pattern_file))
  # calculate Z-scores
  if (z_score == 1) {
    outlist_stats <- calculate_zscores_peakgrouplist(outlist_total)
  }
  # calculate ppm deviation
  outlist_withppm <- calculate_ppmdeviation(outlist_stats)
  #  put columns in correct order
  outlist_ident <- order_columns_peakgrouplist(outlist_withppm)

  # generate output in Excel-readable format:
  remove_columns <- c("mzmin.pgrp", "mzmax.pgrp")
  outlist_ident <- outlist_ident[, -which(colnames(outlist_ident) %in% remove_columns)]
  write.table(outlist_ident, file = paste0("outlist_identified_", scanmode, ".txt"), sep = "\t", row.names = FALSE)
  # generate output in RData format
  save(outlist_ident, file = paste0("outlist_identified_", scanmode, ".RData"))
}
