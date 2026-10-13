# Define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

peakgrouplist_file <- cmd_args[1]
preprocessing_scripts_dir <- cmd_args[2]
thresh <- as.numeric(cmd_args[3])

# Load in function scripts
source(paste0(preprocessing_scripts_dir, "fill_missing_functions.R"))

# Initialize
options(digits = 16)

# Determine scan mode
if (grepl("_pos", peakgrouplist_file)) {
  scanmode <- "positive"
} else if (grepl("_neg", peakgrouplist_file)) {
  scanmode <- "negative"
}

# Get replication pattern for sample names
pattern_file <- paste0(scanmode, "_repl_pattern.RData")
repl_pattern <- get(load(pattern_file))

# Load peak group list and determine output file name
peakgroup_list <- get(load(peakgrouplist_file))

# Replace missing values (zeros) with random noise
peakgrouplist_filled <- fill_missing_intensities(peakgroup_list, repl_pattern, thresh)

# Set name of output file
outputfile_name <- gsub(".RData", "_filled.RData", peakgrouplist_file)

# Save output
save(peakgrouplist_filled, file = outputfile_name)

