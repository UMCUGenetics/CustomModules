#!/usr/bin/Rscript
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
trim_left <- NULL
trim_right <- NULL
breaks_fwhm <- NULL
breaks_fwhm_avg <- NULL
bins <- NULL

# read in mzML file
raw_data <- suppressMessages(xcmsRaw(filepath))

# trim (remove) scans at the start and end
trim_left  <- round(raw_data@scantime[length(raw_data@scantime) * trim])
trim_right <- round(raw_data@scantime[length(raw_data@scantime) * (1 - trim)])

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
save(breaks_fwhm, breaks_fwhm_avg, trim_left, trim_right, file = "breaks.fwhm.RData")
save(high_mz, file = "highest_mz.RData")
