# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

replicate_rdatafile <- cmd_args[1]
resol <- as.numeric(cmd_args[2])
preprocessing_scripts_dir <- cmd_args[3]
# use fixed theshold between noise and signal for peak
peak_thresh <- 2000 

# source functions script
source(paste0(preprocessing_scripts_dir, "peak_finding_functions.R"))

library(dplyr)

# Load output of AssignToBins (peak_list) for a technical replicate
load(replicate_rdatafile)
techrepl_name <- colnames(peak_list$pos)[1]

# load list of technical replicates per sample that passed threshold filter
techreps_passed <- read.table("replicates_per_sample.txt", sep=",")

# Initialize
options(digits = 16)

# run the findPeaks function
scanmodes <- c("positive", "negative")
for (scanmode in scanmodes) {
  # get intensities for scan mode
  if (scanmode == "positive") {
    ints_perscanmode <- peak_list$pos
  } else if (scanmode == "negative") {
    ints_perscanmode <- peak_list$neg
  }
 
  # check whether technical replicate has passed threshold filter for this scanmode
  techreps_scanmode <- techreps_passed[grep(scanmode, techreps_passed[, 3]), ]
  # if techrep is ok, it will be found. If not, skip this techrep.
  if (length(grep(techrepl_name, techreps_scanmode)) == 0) {
    break
  }

  # put mz and intensities into dataframe
  ints_fullrange <- as.data.frame(cbind(mz = as.numeric(rownames(ints_perscanmode)), 
                                        int = as.numeric(ints_perscanmode)))

  # look for m/z range for all peaks
  regions_of_interest <- search_regions_of_interest(ints_fullrange)
 
  # fit Gaussian curve and calculate integrated area under curve
  integrated_peak_df <- integrate_peaks(ints_fullrange, regions_of_interest, resol, peak_thresh)
 
  # add sample name to dataframe
  integrated_peak_df <- as.data.frame(cbind(samplenr = techrepl_name, integrated_peak_df))
 
  # save output to file
  save(integrated_peak_df, file = paste0(techrepl_name, "_", scanmode, ".RData"))
}

