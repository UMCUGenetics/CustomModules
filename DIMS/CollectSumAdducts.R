#!/usr/bin/Rscript
## adapted from 12-collectSamplesAdded.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

# outdir <- cmd_args[1]
# scanmode <- cmd_args[2]
scanmode <- c("positive", "negative")

for (scanmode in scanmodes) {
  # collect all AdductSums part files for each scanmode
  adductsum_part_files <- list.files("./", pattern = scanmode)

  outlist.tot <- NULL
  for (i in 1:length(adductsum_part_files)) {
    load(adductsum_part_files[i])
    outlist.tot <- rbind(outlist.tot, adductsum)
  }

  # save output file
  save(outlist.tot, file="/AdductSums_", scanmode, ".RData")
}

