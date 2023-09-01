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
  filled_files <- list.files("./", full.names=TRUE, pattern=scanmode)
  # load files and combine into one object
  outlist.tot <- NULL
  for (i in 1:length(filled_files)) {
    load(filled_files[i])
    print(filled_files[i])
    outlist.tot <- rbind(outlist.tot, peakgrouplist_filled)
  }

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
    nr.removed.samples <- length(which(repl_pattern[] == "character(0)"))
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

  # filter identified compounds
  index.1 <- which((outlist.tot[,"assi_HMDB"]!="") & (!is.na(outlist.tot[,"assi_HMDB"])))
  index.2 <- which((outlist.tot[,"iso_HMDB"]!="") & (!is.na(outlist.tot[,"iso_HMDB"])))
  index <- union(index.1,index.2)
  outlist.ident <- outlist.tot[index,]
  outlist.not.ident <- outlist.tot[-index,]

  if (z_score == 1) {
    outlist.ident$ppmdev <- as.numeric(outlist.ident$ppmdev)
    outlist.ident <- outlist.ident[which(outlist.ident["ppmdev"] >= -ppm & outlist.ident["ppmdev"] <= ppm),]
  }
  # NAs in theormz_noise 
  outlist.ident$theormz_noise[which(is.na(outlist.ident$theormz_noise))] <- 0
  outlist.ident$theormz_noise <- as.numeric(outlist.ident$theormz_noise)
  outlist.ident$theormz_noise[which(is.na(outlist.ident$theormz_noise))] <- 0
  outlist.ident$theormz_noise <- as.numeric(outlist.ident$theormz_noise)

  # save output for identified peak groups (not.identified later)
  # save(outlist.not.ident, file="./outlist_not_identified_", scanmode, ".RData")
  save(outlist.ident, file=paste0("./outlist_identified_", scanmode, ".RData"))

}
