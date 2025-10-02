library(dplyr)

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

sample_name <- cmd_args[1]
techreps <- cmd_args[2]
scanmode <- cmd_args[3]
preprocessing_scripts_dir <- cmd_args[4]
tech_reps <- strsplit(techreps, ";")[[1]]

# load in function scripts
source(paste0(preprocessing_scripts_dir, "average_peaks_functions.R"))

# set ppm as fixed value, not the same ppm as in peak grouping
ppm_peak <- 2

# Initialize per sample
peaklist_allrepl <- NULL
nr_repl_persample <- 0
averaged_peaks <- matrix(0, nrow = 0, ncol = 6) 
colnames(averaged_peaks) <- c("samplenr", "mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")

# load RData files of technical replicates belonging to biological sample
for (file_nr in 1:length(tech_reps)) {
  tech_repl_file <- paste0(tech_reps[file_nr], "_", scanmode, ".RData")
  tech_repl <- get(load(tech_repl_file))
  # combine data for all technical replicates
  peaklist_allrepl <- rbind(peaklist_allrepl, tech_repl)
}
# sort on mass
peaklist_allrepl_df <- as.data.frame(peaklist_allrepl)
peaklist_allrepl_df$mzmed.pkt <- as.numeric(peaklist_allrepl_df$mzmed.pkt) 
peaklist_allrepl_df$height.pkt <- as.numeric(peaklist_allrepl_df$height.pkt) 
peaklist_allrepl_sorted <- peaklist_allrepl_df %>% arrange(mzmed.pkt)

# average over technical replicates
averaged_peaks <- average_peaks_per_sample(peaklist_allrepl_sorted)
save(averaged_peaks, file = paste0("AvgPeaks_", sample_name, "_", scanmode, ".RData"))

