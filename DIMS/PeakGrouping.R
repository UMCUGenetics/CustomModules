# load required packages
library("dplyr")
library("argparse")

parser <- ArgumentParser(description = "PeakGrouping")

parser$add_argument("--hmdbpart_file", dest = "hmdbpart_file",
                    help = "RData file containing part of HMDB database (plus adducts and isotopes)", required = TRUE)
parser$add_argument("--preprocessing_scripts_dir", dest = "preprocessing_scripts_dir",
                    help = "File path to the directory containing functions used", required = TRUE)
parser$add_argument("--ppm", dest = "ppm",
                    help = "Value for ppm deviation allowed in peak group annotation", required = TRUE)

args <- parser$parse_args()

hmdbpart_file <- args$hmdbpart_file
ppm <- as.numeric(args$ppm)

# load in function scripts
source(paste0(args$preprocessing_scripts_dir, "peak_grouping_functions.R"))

options(digits = 16)

# load part of the HMDB
hmdb_add_iso <- get(load(args$hmdbpart_file))
# make sure any internal standard entries are at the top of the dataframe
internal_standard_indices <- grep("\\(IS\\)", hmdb_add_iso$CompoundName)
if (length(internal_standard_indices) > 0) {
  hmdb_add_iso <- rbind(hmdb_add_iso[internal_standard_indices, ], hmdb_add_iso[-internal_standard_indices, ])
}

# determine appropriate scanmode based on hmdbpart_file
if (grepl("negative", basename(hmdbpart_file))) {
  scanmode <- "negative"
  column_label <- "MNeg"
} else if (grepl("positive", basename(hmdbpart_file))) {
  scanmode <- "positive"
  column_label <- "Mpos"
}

# determine batch number of HMDB part file
batch_number <- strsplit(basename(hmdbpart_file), ".", fixed = TRUE)[[1]][2]

# load file with averaged peaks
avg_peaks_file <- paste0("AvgPeaks_", scanmode, ".RData")
load(avg_peaks_file)
outlist_df <- as.data.frame(outlist_total)
outlist_df$mzmed.pkt <- as.numeric(outlist_df$mzmed.pkt)
outlist_df$height.pkt <- as.numeric(outlist_df$height.pkt)
rm(outlist_total)
sample_names <- unique(outlist_df$samplenr)

## peak grouping
peakgrouplist <- NULL
# limit the peaklist to the m/z range in the HMDB part, with ppm tolerance
minmz_hmdbpart <- min(hmdb_add_iso[, column_label])
maxmz_hmdbpart <- max(hmdb_add_iso[, column_label])
mz_tolerance <- (maxmz_hmdbpart * ppm) / 10^6
outlist_mzrange <- outlist_df[outlist_df$mzmed.pkt > (minmz_hmdbpart - mz_tolerance) &
                              outlist_df$mzmed.pkt < (maxmz_hmdbpart + mz_tolerance), ]
# sort by descending intensity
outlist_sorted <- outlist_mzrange %>% dplyr::arrange(desc(height.pkt))

# find peak groups
ints_sorted <- find_peak_groups(outlist_sorted, mz_tolerance, sample_names)

# do annotation
peakgrouplist_identified <- annotate_peak_groups(ints_sorted, hmdb_add_iso, column_label, mz_tolerance)

# write output to file
save(peakgrouplist_identified, file = paste0(batch_number, "_", scanmode, "_identified.RData"))

