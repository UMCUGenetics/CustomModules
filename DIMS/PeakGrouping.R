# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

# load required packages
library("dplyr")

hmdb_part_file <- cmd_args[1]
preprocessing_scripts_dir <- cmd_args[2]
ppm <- as.numeric(cmd_args[3])

# load in function scripts
source(paste0(preprocessing_scripts_dir, "peak_grouping_functions.R"))

options(digits = 16)

# load part of the HMDB
hmdb_add_iso <- get(load(hmdb_part_file))

# determine appropriate scanmode based on hmdb_part_file
if (grepl("negative", basename(hmdb_part_file))) {
  scanmode <- "negative"
  column_label <- "MNeg"
} else if (grepl("positive", basename(hmdb_part_file))) {
  scanmode <- "positive"
  column_label <- "Mpos"
}

# determine batch number of HMDB part file
batch_number <- strsplit(basename(hmdb_part_file), ".", fixed = TRUE)[[1]][2]

# load file with spectrum peaks
spec_peaks_file <- paste0("SpectrumPeaks_", scanmode, ".RData")
load(spec_peaks_file)
outlist_df <- as.data.frame(outlist_total)
rm(outlist_total)
sample_names <- unique(outlist_df$samplenr)

## peak grouping
peakgrouplist <- NULL
# limit the peaklist to the m/z range in the HMDB part, with ppm tolerance
minmz_hmdbpart <- min(hmdb_add_iso[ , column_label])
maxmz_hmdbpart <- max(hmdb_add_iso[ , column_label])
mz_tolerance <- (maxmz_hmdbpart * ppm) / 10^6
outlist_mzrange <- outlist_df[outlist_df$mzmed.pkt > (minmz_hmdbpart - mz_tolerance) &
                              outlist_df$mzmed.pkt < (maxmz_hmdbpart + mz_tolerance), ]
# sort by descending intensity
outlist_mzrange$mzmed.pkt <- as.numeric(outlist_mzrange$mzmed.pkt)  # lapply?
outlist_mzrange$height.pkt <- as.numeric(outlist_mzrange$height.pkt)
outlist_sorted <- outlist_mzrange %>% dplyr::arrange(desc(height.pkt))

# find peak groups
ints_sorted <- find_peak_groups(outlist_sorted, mz_tolerance, sample_names)

# count the number of non-zero intensities per row. First 3 columns are m/z
nrsamples <- apply(ints_sorted[, 4:ncol(ints_sorted)], 1, function(x) sum(x > 0))

# do annotation
peakgrouplist_identified <- annotate_peak_groups(ints_sorted, hmdb_add_iso)

# write output to file
save(peakgrouplist_identified, file = paste0(batch_number, "_", scanmode, "_identified.RData"))

