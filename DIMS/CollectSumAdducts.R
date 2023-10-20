#!/usr/bin/Rscript
## adapted from 12-collectSamplesAdded.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)

scanmodes <- c("positive", "negative")

for (scanmode in scanmodes) {
  # collect all AdductSums part files for each scanmode
  adductsum_part_files <- list.files("./", pattern = scanmode)

  outlist.tot <- NULL
  for (i in 1:length(adductsum_part_files)) {
    load(adductsum_part_files[i])
    outlist.tot <- rbind(outlist.tot, adductsum)
  }

  # save output file
  save(outlist.tot, file=paste0("AdductSums_", scanmode, ".RData"))
}

