# define parameters
# ppm as fixed value, not the same ppm as in peak grouping
ppm_peak <- 2

library(dplyr)

scanmodes <- c("positive", "negative")

for (scanmode in scanmodes){
  # get sample names
  load(paste0(scanmode, "_repl_pattern.RData"))
  sample_names <- names(repl_pattern_filtered)
  # initialize
  outlist_total <- NULL
  # for each biological sample, average peaks in technical replicates
  for (sample_name in sample_names) {
    print(sample_name)
    # Initialize per sample
    peaklist_allrepl <- NULL
    nr_repl_persample <- 0
    # averaged_peaks <- matrix(0, nrow = 10 ^ 7, ncol = 6) # how big does it need to be?
    averaged_peaks <- matrix(0, nrow = 0, ncol = 6) # append
    colnames(averaged_peaks) <- c("samplenr", "mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")
    # load RData files of technical replicates belonging to biological sample
    sample_techreps_file <- repl_pattern_filtered[sample_name][[1]]
    for (file_nr in 1:length(sample_techreps_file)) {
      print(sample_techreps_file[file_nr])
      tech_repl_file <- paste0(sample_techreps_file[file_nr], "_positive.RData")
      tech_repl <- get(load(tech_repl_file))
      # combine data for all technical replicates
      peaklist_allrepl <- rbind(peaklist_allrepl, tech_repl)
      # count number of replicates for each biological sample
      nr_repl_persample <- nr_repl_persample + 1
    }
    # sort on mass
    peaklist_allrepl_df <- as.data.frame(peaklist_allrepl)
    peaklist_allrepl_df$mzmed.pkt <- as.numeric(peaklist_allrepl_df$mzmed.pkt) 
    peaklist_allrepl_df$height.pkt <- as.numeric(peaklist_allrepl_df$height.pkt) 
    # peaklist_allrepl_sorted <- peaklist_allrepl_df %>% arrange(desc(height.pkt))
    peaklist_allrepl_sorted <- peaklist_allrepl_df %>% arrange(mzmed.pkt)
    # average over technical replicates
    while (nrow(peaklist_allrepl_sorted) > 1) {
      # store row numbers
      peaklist_allrepl_sorted$rownr <- 1:nrow(peaklist_allrepl_sorted)
      # find the peaks in the dataset with corresponding m/z plus or minus tolerance
      reference_mass <- peaklist_allrepl_sorted$mzmed.pkt[1]
      mz_tolerance <- (reference_mass * ppm_peak) / 10^6
      minmz_ref <- reference_mass - mz_tolerance
      maxmz_ref <- reference_mass + mz_tolerance
      select_peak_indices <- which((peaklist_allrepl_sorted$mzmed.pkt > minmz_ref) & (peaklist_allrepl_sorted$mzmed.pkt < maxmz_ref))
      select_peaks <- peaklist_allrepl_sorted[select_peak_indices, ]
      nrsamples <- length(select_peak_indices)
      # put averaged intensities into a new row and append to averaged_peaks
      averaged_1peak <- matrix(0, nrow = 1, ncol = 6) 
      colnames(averaged_1peak) <- c("samplenr", "mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")
      # calculate m/z values for peak group
      averaged_1peak[1, "mzmed.pkt"] <- mean(select_peaks$mzmed.pkt)
      averaged_1peak[1, "mzmin.pkt"] <- min(select_peaks$mzmed.pkt)
      averaged_1peak[1, "mzmax.pkt"] <- max(select_peaks$mzmed.pkt)
      averaged_1peak[1, "fq"] <- nrsamples
      averaged_1peak[1, "height.pkt"] <- mean(select_peaks$height.pkt)
      # put intensities into proper columns
      peaklist_allrepl_sorted <- peaklist_allrepl_sorted[-select_peaks$rownr, ]
      averaged_peaks <- rbind(averaged_peaks, averaged_1peak)
    }
    # add sample name to first column and append to outlist_total for all samples
    averaged_peaks[ , "samplenr"] <- sample_name
    outlist_total <- rbind(outlist_total, averaged_peaks)
  }
  save(outlist_total, file = paste0("AvgPeaks_", scanmode, ".RData"))
}

