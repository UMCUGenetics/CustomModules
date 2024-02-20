#!/usr/bin/Rscript
## adapted from 10-collectSamplesFilled.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

scripts_dir        <- cmd_args[1]
ppm     <- as.numeric(cmd_args[2])
z_score <- as.numeric(cmd_args[3])

source(paste0(scripts_dir, "AddOnFunctions/mergeDuplicatedRows.R"))
source(paste0(scripts_dir, "AddOnFunctions/statistics_z.R"))

# for each scan mode, collect all filled peak group lists
scanmodes <- c("positive", "negative")

for (scanmode in scanmodes) {
  # get list of files
  filled_file <- paste0("./PeakGroupList_", scanmode, "_Unidentified_filled.RData")
  # load file
  outlist_total <- get(load(filled_file))

  # remove duplicates; peak groups with exactly the same m/z
  outlist_total <- mergeDuplicatedRows(outlist_total)

  # sort on mass
  outlist_total <- outlist_total[order(outlist_total[, "mzmed.pgrp"]), ]

  # load replication pattern
  pattern_file <- paste0(scanmode, "_repl_pattern.RData")
  repl_pattern <- get(load(pattern_file))

  if (z_score == 1) {
    # calculate Z-scores
    outlist_stats <- statistics_z(outlist_total, sortCol = NULL, adducts = FALSE)
    nr_removed_samples <- length(which(repl_pattern[] == "character(0)"))
    order_index_int <- order(colnames(outlist_stats)[8:(length(repl_pattern) - nr_removed_samples + 7)])
    outlist_stats_more <- cbind(outlist_stats[, 1:7],
                                outlist_stats[, (length(repl_pattern) - nr_removed_samples + 8):(length(repl_pattern) - nr_removed_samples + 8 + 6)],
                                outlist_stats[, 8:(length(repl_pattern) - nr_removed_samples + 7)][order_index_int],
                                outlist_stats[, (length(repl_pattern) - nr_removed_samples + 5 + 10):ncol(outlist_stats)])

    tmp_index <- grep("_Zscore", colnames(outlist_stats_more), fixed = TRUE)
    tmp_index_order <- order(colnames(outlist_stats_more[, tmp_index]))
    tmp <- outlist_stats_more[, tmp_index[tmp_index_order]]
    outlist_stats_more <- outlist_stats_more[, -tmp_index]
    outlist_stats_more <- cbind(outlist_stats_more, tmp)
    outlist_total <- outlist_stats_more
  }

  outlist_not_ident <- outlist_total

  # Extra output in Excel-readable format:
  remove_columns <- c("fq.best", "fq.worst", "mzmin.pgrp", "mzmax.pgrp")
  remove_colindex <- which(colnames(outlist_not_ident) %in% remove_columns)
  outlist_not_ident <- outlist_not_ident[, -remove_colindex]
  write.table(outlist_not_ident, file = paste0("unidentified_outlist_", scanmode, ".txt"), sep = "\t", row.names = FALSE)

}
