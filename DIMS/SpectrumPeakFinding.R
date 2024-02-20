#!/usr/bin/Rscript
# adapted from 5-collectSamples.R

# define parameters
scanmodes <- c("positive", "negative")

# Check whether all jobs terminated correct!
not_run <- NULL

# collect spectrum peaks for each scanmode
for (scanmode in scanmodes) {
  # load peak lists of all biological samples
  peaklist_files <- list.files("/.", full.names = TRUE, pattern = paste0("*_", scanmode, ".RData"))

  # get sample names
  load(paste0("./", scanmode, "_repl_pattern", ".RData"))
  group_names <- names(repl_pattern_filtered)
  for (sample_nr in 1:length(group_names)) {
    group <- paste0(input_dir, "/", paste0(paste(group_names[sample_nr], scanmode, sep = "_"), ".RData"))
    if (!(group %in% peaklist_files)) {
      not_run <- c(not_run, group)
    }
  }

  # Collecting samples
  outlist_total <- NULL
  for (file_nr in 1:length(peaklist_files)) {
    cat("\n", peaklist_files[file_nr])
    load(peaklist_files[file_nr])
    if (is.null(outlist.persample) || (dim(outlist.persample)[1] == 0)) {
      tmp <- strsplit(peaklist_files[file_nr], "/")[[1]]
      fname <- tmp[length(tmp)]
      fname <- strsplit(fname, ".RData")[[1]][1]
      fname <- substr(fname, 13, nchar(fname))
      if (file_nr == 1) {
        outlist_total <- c(fname, rep("-1", 5))
      } else {
        outlist_total <- rbind(outlist_total, c(fname, rep("-1", 5)))
      }
    } else {
      if (file_nr == 1) {
        outlist_total <- outlist.persample
      } else {
        outlist_total <- rbind(outlist_total, outlist.persample)
      }
    }
  }

  # remove negative values
  index <- which(outlist_total[, "height.pkt"] <= 0)
  if (length(index) > 0) outlist_total <- outlist_total[-index, ]
  index <- which(outlist_total[, "mzmed.pkt"] <= 0)
  if (length(index) > 0) outlist_total <- outlist_total[-index, ]

  save(outlist_total, file = paste0("./SpectrumPeaks_", scanmode, ".RData"))

  if (!is.null(not_run)) {
    for (i in 1:length(not_run)) {
      message(paste(not_run[i], "was not generated"))
    }
  }
}
