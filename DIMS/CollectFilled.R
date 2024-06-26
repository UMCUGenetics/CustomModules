## adapted from 10-collectSamplesFilled.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

scripts_dir <- cmd_args[1]
ppm <- as.numeric(cmd_args[2])
z_score <- as.numeric(cmd_args[3])

source(paste0(scripts_dir, "merge_duplicate_rows.R"))
source(paste0(scripts_dir, "calculate_zscores.R"))

# for each scan mode, collect all filled peak group lists
scanmodes <- c("positive", "negative")
for (scanmode in scanmodes) {
  # get list of files
  filled_files <- list.files("./", full.names = TRUE, pattern = paste0(scanmode, "_identified_filled"))
  # load files and combine into one object
  outlist_total <- NULL
  for (file_nr in 1:length(filled_files)) {
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
    outlist_stats <- calculate_zscores(outlist_total, adducts = FALSE)
    nr_removed_samples <- length(which(repl_pattern[] == "character(0)"))
    order_index_int <- order(colnames(outlist_stats)[8:(length(repl_pattern) - nr_removed_samples + 7)])
    outlist_stats_more <- cbind(
      outlist_stats[, 1:7],
      outlist_stats[, (length(repl_pattern) - nr_removed_samples + 8):(length(repl_pattern) - nr_removed_samples + 8 + 6)],
      outlist_stats[, 8:(length(repl_pattern) - nr_removed_samples + 7)][order_index_int],
      outlist_stats[, (length(repl_pattern) - nr_removed_samples + 5 + 10):ncol(outlist_stats)]
    )
    # sort Z-score columns and append to peak group list
    tmp_index <- grep("_Zscore", colnames(outlist_stats_more), fixed = TRUE)
    tmp_index_order <- order(colnames(outlist_stats_more[, tmp_index]))
    tmp <- outlist_stats_more[, tmp_index[tmp_index_order]]
    outlist_stats_more <- outlist_stats_more[, -tmp_index]
    outlist_stats_more <- cbind(outlist_stats_more, tmp)
    outlist_total <- outlist_stats_more
  }
  
  # make a copy of the outlist
  outlist_ident <- outlist_total
  # select identified peak groups if ppm deviation is within limits
  if (z_score == 1) {
    outlist_ident$ppmdev <- as.numeric(outlist_ident$ppmdev)
    outlist_ident <- outlist_ident[which(outlist_ident["ppmdev"] >= -ppm & outlist_ident["ppmdev"] <= ppm), ]
  }
  # take care of NAs in theormz_noise
  outlist_ident$theormz_noise[which(is.na(outlist_ident$theormz_noise))] <- 0
  outlist_ident$theormz_noise <- as.numeric(outlist_ident$theormz_noise)
  outlist_ident$theormz_noise[which(is.na(outlist_ident$theormz_noise))] <- 0
  outlist_ident$theormz_noise <- as.numeric(outlist_ident$theormz_noise)

  # Extra output in Excel-readable format:
  remove_columns <- c("fq.best", "fq.worst", "mzmin.pgrp", "mzmax.pgrp")
  remove_colindex <- which(colnames(outlist_ident) %in% remove_columns)
  outlist_ident <- outlist_ident[, -remove_colindex]
  write.table(outlist_ident, file = paste0("outlist_identified_", scanmode, ".txt"), sep = "\t", row.names = FALSE)
  # output in RData format
  save(outlist_ident, file = paste0("outlist_identified_", scanmode, ".RData"))
}
