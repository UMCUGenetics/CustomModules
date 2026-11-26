# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

scripts_dir <- cmd_args[1]

# for each scan mode, collect all averaged peak lists per biological sample
scanmodes <- c("positive", "negative")
for (scanmode in scanmodes) {
  # get list of files
  filled_files <- list.files("./", full.names = TRUE, pattern = paste0(scanmode, ".RData"))
  # load files and combine into one object
  outlist_total <- NULL
  for (file_nr in 1:length(filled_files)) {
    peaklist_averaged <- get(load(filled_files[file_nr]))
    outlist_total <- rbind(outlist_total, peaklist_averaged)
  }
  save(outlist_total, file = paste0("AvgPeaks_", scanmode, ".RData"))
}
