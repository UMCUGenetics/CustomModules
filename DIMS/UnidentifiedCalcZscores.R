#!/usr/bin/Rscript
## adapted from 10-collectSamplesFilled.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n")

scripts_dir        <- cmd_args[1]
ppm     <- as.numeric(cmd_args[2])
z_score <- as.numeric(cmd_args[3])

source(paste0(scripts_dir, "AddOnFunctions/mergeDuplicatedRows.R"))
source(paste0(scripts_dir, "AddOnFunctions/statistics_z.R"))
# source(paste0(scripts_dir, "AddOnFunctions/normalization_2.1.R"))

# for each scan mode, collect all filled peak group lists
scanmodes <- c("positive", "negative")

for (scanmode in scanmodes) {
  # get list of files
  # filled_files <- list.files("./", full.names=TRUE, pattern=scanmode)
  filled_file <- paste0("./PeakGroupList_", scanmode, "_Unidentified_filled.RData")
  print(filled_file)
  # load file 
  outlist.tot <- get(load(filled_file))

  # remove duplicates; peak groups with exactly the same m/z
  outlist.tot <- mergeDuplicatedRows(outlist.tot)

  # sort on mass
  outlist.tot <- outlist.tot[order(outlist.tot[ ,"mzmed.pgrp"]),]

  # load replication pattern
  pattern_file <- paste0(scanmode, "_repl_pattern.RData")
  repl_pattern <- get(load(pattern_file))

  # Normalization: not done.
  # if (normalization != "disabled") {
  #   outlist.tot = normalization_2.1(outlist.tot, fileName, names(repl.pattern.filtered), on=normalization, assi_label="assi_HMDB")
  # }

  if (z_score == 1) {
    outlist.stats <- statistics_z(outlist.tot, sortCol=NULL, adducts=FALSE)
    nr.removed.samples <- length(which(repl_pattern[]=="character(0)"))
    order.index.int <- order(colnames(outlist.stats)[8:(length(repl_pattern)-nr.removed.samples+7)])
    outlist.stats.more <- cbind(outlist.stats[,1:7],
                               outlist.stats[,(length(repl_pattern)-nr.removed.samples+8):(length(repl_pattern)-nr.removed.samples+8+6)],
                               outlist.stats[,8:(length(repl_pattern)-nr.removed.samples+7)][order.index.int],
                               outlist.stats[,(length(repl_pattern)-nr.removed.samples+5+10):ncol(outlist.stats)])
  
    tmp.index <- grep("_Zscore", colnames(outlist.stats.more), fixed = TRUE)
    tmp.index.order <- order(colnames(outlist.stats.more[,tmp.index]))
    tmp <- outlist.stats.more[,tmp.index[tmp.index.order]]
    outlist.stats.more <- outlist.stats.more[,-tmp.index]
    outlist.stats.more <- cbind(outlist.stats.more,tmp)
    outlist.tot <- outlist.stats.more
  }

  outlist.not.ident = outlist.tot

  # Extra output in Excel-readable format:
  remove_columns <- c("fq.best", "fq.worst", "mzmin.pgrp", "mzmax.pgrp")
  remove_colindex <- which(colnames(outlist.not.ident) %in% remove_columns)
  outlist.not.ident <- outlist.not.ident[ , -remove_colindex]
  write.table(outlist.not.ident, file=paste0("unidentified_outlist_", scanmode, ".txt"), sep="\t", row.names = FALSE)

}
