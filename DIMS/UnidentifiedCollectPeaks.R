#!/usr/bin/Rscript
## adapted from 7-collectSamplesGroupedHMDB.R

options(digits = 16)

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

ppm <- as.numeric(cmd_args[1])

scanmodes <- c("positive", "negative")

for (scanmode in scanmodes) {
  # get list of all files that contain lists of peaks that were used in identified peak grouping
  files <- list.files("./", pattern = paste(scanmode, "_peaks_used.RData", sep = ""))
  # load the list of all peaks
  load(paste0("SpectrumPeaks_", scanmode, ".RData"))

  # Make a list of indexes of peaks that have been identified, then remove these from the peaklist.
  remove <- NULL
  for (file_index in 1:length(files)) {
    # load list_of_peaks_used_in_peak_groups_identified
    load(files[file_index])
    remove <- c(remove, 
		which(outlist_total[, "mzmed.pkt"] %in% list_of_peaks_used_in_peak_groups_identified[, "mzmed.pkt"])) # nolint
  }
  outlist_rest <- outlist_total[-remove, ]

  # sort on mass
  outlist <- outlist_rest[order(as.numeric(outlist_rest[, "mzmed.pkt"])), ]
  # save output
  # save(outlist, file = paste0("SpectrumPeaks_", scanmode, "_Unidentified.RData"))

  size_parts <- 10000
  # start <- 1

  num_parts <- ceiling(nrow(outlist) / size_parts)

  for (part in 1:num_parts){
    if (part == 1) {
      start_part <- 1
      end_part <- size_parts
    } else {
      start_part <- (part - 1) * size_parts + 1
    }
    if (part == num_parts) {
      end_part <- nrow(outlist)
    } else if (part != 1) {
      end_part <- part * size_parts
    }

    outlist_part <- outlist[c(start_part:end_part), ]
    # # add ppm extra before start
    # if (part != 1) {
    #   mz_start <- outlist_part[1, "mzmed.pkt"]
    #   mz_ppm_range <- ppm * as.numeric(mz_start) / 1e+06
    #   mz_start_min_ppm <- mz_start - mz_ppm_range
    #   outlist_before_part <- outlist %>% filter(mzmed.pkt >= mz_start_min_ppm & mzmed.pkt < mz_start)

    #   outlist_part <- rbind(outlist_before_part, outlist_part)
    # }

    # add ppm extra after end
    if (part != num_parts) {
      mz_end <- outlist_part[nrow(outlist_part), "mzmed.pkt"]
      mz_ppm_range <- ppm * as.numeric(mz_end) / 1e+06
      mz_end_plus_ppm <- mz_end + mz_ppm_range
      outlist_after_part <- outlist %>% filter(mzmed.pkt > mz_end & mzmed.pkt <= mz_end_plus_ppm)

      outlist_part <- rbind(outlist_part, outlist_after_part)
    }

    save(outlist_part, file = paste0(scanmode, "_", paste0("unidentified_part_", part, ".RData")))
  }
}
