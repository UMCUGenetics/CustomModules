#!/usr/bin/Rscript
## adapted from 7-collectSamplesGroupedHMDB.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

ppm <- as.numeric(cmd_args[1])
outdir <- "./"

scanmodes <- c("positive", "negative")

for (scanmode in scanmodes) {
  # get list of all files that contain lists of peaks that were used in identified peak grouping
  files <- list.files("./", pattern = paste(scanmode, "_peaks_used.RData", sep = ""))
  # load the list of all peaks
  load(paste0("SpectrumPeaks_", scanmode, ".RData"))

  # Make a list of indexes of peaks that have been identified, then remove these from the peaklist.
  remove <- NULL
  for (i in 1:length(files)) {
    # load list_of_peaks_used_in_peak_groups_identified
    load(files[i])
    remove <- c(remove, which(outlist_total[, "mzmed.pkt"] %in% list_of_peaks_used_in_peak_groups_identified[i, "mzmed.pkt"]))
  }
  outlist_rest <- outlist_total[-remove, ]

  # sort on mass
  outlist <- outlist_rest[order(as.numeric(outlist_rest[, "mzmed.pkt"])), ]
  # save output
  save(outlist, file = paste0(outdir, "/SpectrumPeaks_", scanmode, "_Unidentified.RData"))

}
