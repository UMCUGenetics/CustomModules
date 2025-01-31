## adapted from 1-generateBreaksFwhm.HPC.R ##

# load required package
suppressPackageStartupMessages(library("xcms"))

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

filepath <- cmd_args[1]
outdir <- cmd_args[2]
trim <- as.numeric(cmd_args[3])
resol <- as.numeric(cmd_args[4])

# initialize
trim_left_pos <- NULL
trim_right_pos <- NULL
trim_left_neg <- NULL
trim_right_neg <- NULL
breaks_fwhm <- NULL
breaks_fwhm_avg <- NULL
bins <- NULL

# read in mzML file
raw_data <- suppressMessages(xcms::xcmsRaw(filepath))

# Get time values for positive and negative scans
pos_times <- raw_data@scantime[raw_data@polarity == "positive"]
neg_times <- raw_data@scantime[raw_data@polarity == "negative"]

# trim (remove) scans at the start and end for positive
trim_left_pos  <- round(pos_times[length(pos_times) * (trim * 1.5)]) # 15% aan het begin
trim_right_pos <- round(pos_times[length(pos_times) * (1 - (trim * 0.5))]) # 5% aan het eind

# trim (remove) scans at the start and end for negative
trim_left_neg  <- round(neg_times[length(neg_times) * trim])
trim_right_neg <- round(neg_times[length(neg_times) * (1 - trim)])

# Mass range m/z
low_mz  <- raw_data@mzrange[1]
high_mz <- raw_data@mzrange[2]

# determine number of segments (bins)
nr_segments <- 2 * (high_mz - low_mz)
segment <- seq(from = low_mz, to = high_mz, length.out = nr_segments + 1)

# determine start and end of each bin.
for (i in 1:nr_segments) {
  start_segment <- segment[i]
  end_segment <- segment[i+1]
  resol_mz <- resol * (1 / sqrt(2) ^ (log2(start_segment / 200)))
  fwhm_segment <- start_segment / resol_mz
  breaks_fwhm <- c(breaks_fwhm, seq(from = (start_segment + fwhm_segment), to = end_segment, by = 0.2 * fwhm_segment))
  # average the m/z instead of start value
  range <- seq(from = (start_segment + fwhm_segment), to = end_segment, by = 0.2 * fwhm_segment)
  delta_mz <- range[2] - range[1]
  breaks_fwhm_avg <- c(breaks_fwhm_avg, range + 0.5 * delta_mz)
}

# generate output file
save(breaks_fwhm, breaks_fwhm_avg, file = "breaks.fwhm.RData")
save(trim_left_pos, trim_right_pos, trim_left_neg, trim_right_neg, file = "trim_params.RData")
save(high_mz, file = "highest_mz.RData")
