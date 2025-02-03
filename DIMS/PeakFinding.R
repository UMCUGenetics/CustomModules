# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

replicate_mzmlfile <- cmd_args[1]
breaks_file <- cmd_args[2]
resol <- as.numeric(cmd_args[3])
scripts_dir <- cmd_args[4]
peak_thresh <- 2000

# load in function scripts
preprocessing_scripts_dir <- gsub("Utils", "preprocessing", scripts_dir)
source(paste0(preprocessing_scripts_dir, "peak_finding_functions.R"))

load(breaks_file)

# Load output of AssignToBins for a sample
sample_techrepl <- get(load(replicate_mzmlfile))
scanmodes <- c("positive", "negative")
ints_pos <- peak_list$pos
ints_neg <- peak_list$neg

# Initialize
options(digits = 16)

# run the findPeaks function
techrepl_name <- colnames(ints_pos)[1]

for (scanmode in scanmodes) {
  # turn dataframe with intensities into a named list
  if (scanmode == "positive") {
    ints_perscanmode <- peak_list$pos
  } else if (scanmode == "negative") {
    ints_perscanmode <- peak_list$neg
  }
  
  ints_fullrange <- as.vector(ints_perscanmode)
  names(ints_fullrange) <- rownames(ints_perscanmode)
  
  start.time <- Sys.time()
  
  # look for m/z range for all peaks
  allpeaks_values <- search_mzrange(ints_fullrange, resol, techrepl_name, scanmode, peak_thresh)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  # turn the list into a dataframe
  outlist_persample <- NULL
  outlist_persample <- cbind("samplenr" = allpeaks_values$nr,
                             "mzmed.pkt" = allpeaks_values$mean,
                             "fq" = allpeaks_values$qual,
                             "mzmin.pkt" = allpeaks_values$min,
                             "mzmax.pkt" = allpeaks_values$max,
                             "height.pkt" = allpeaks_values$area)
  
  # remove peaks with height = 0
  outlist_persample <- outlist_persample[outlist_persample[, "height.pkt"] != 0, ]
  
  # save output to file
  save(outlist_persample, file = paste0(techrepl_name, "_", scanmode, ".RData"))
  
  # generate text output to log file on number of spikes for this sample
  # spikes are peaks that are too narrow, e.g. 1 data point
  cat(paste("There were", allpeaks_values$spikes, "spikes"))
  
}

