# define parameters
library("argparse")

parser <- ArgumentParser(description = "FillMissing")

parser$add_argument("--peakgrouplist_file", dest = "peakgrouplist_file",
                    help = "Peak group list file", required = TRUE)
parser$add_argument("--preprocessing_scripts_dir", dest = "preprocessing_scripts_dir",
                    help = "File path to the directory containing functions used", required = TRUE)
parser$add_argument("--thresh_noise", dest = "thresh_noise",
                    help = "Threshold value for noise on peak level (typically 2000)", required = TRUE)

args <- parser$parse_args()

peakgrouplist_file <- args$peakgrouplist_file
thresh <- as.numeric(args$thresh_noise)

# load in function scripts
source(paste0(args$preprocessing_scripts_dir, "fill_missing_functions.R"))

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
