# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

peakgrouplist_file <- cmd_args[1]
preprocessing_scripts_dir <- cmd_args[2]
thresh <- as.numeric(cmd_args[3])
resol <- as.numeric(cmd_args[4])
ppm <- as.numeric(cmd_args[5])

# load in function scripts
source(paste0(preprocessing_scripts_dir, "fill_missing_functions.R"))

# determine scan mode
if (grepl("_pos", peakgrouplist_file)) {
  scanmode <- "positive"
} else if (grepl("_neg", peakgrouplist_file)) {
  scanmode <- "negative"
}

# get replication pattern for sample names
pattern_file <- paste0(scanmode, "_repl_pattern.RData")
repl_pattern <- get(load(pattern_file))

# load peak group list and determine output file name
peakgroup_list <- get(load(peakgrouplist_file))

# replace missing values (zeros) with random noise
peakgrouplist_filled <- fill_missing_intensities(peakgroup_list, repl_pattern, thresh)

# set name of output file
outputfile_name <- gsub(".RData", "_filled.RData", peakgrouplist_file)

# save output
save(peakgrouplist_filled, file = outputfile_name)
