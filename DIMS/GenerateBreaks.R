# load required package
suppressPackageStartupMessages(library("xcms"))

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

filepath <- cmd_args[1]
trim <- as.numeric(cmd_args[2])
resol <- as.numeric(cmd_args[3])

# read in mzML file
raw_data <- suppressMessages(xcms::xcmsRaw(filepath))

# get trim parameters and save them to file
get_trim_parameters(raw_data@scantime, raw_data@polarity)

# create breaks of bins for intensities. Bin size is a function of fwhm which is a function of m/z
get_breaks_for_bins(raw_data$mzrange, resol)

# Determine maximum m/z and save to file
high_mz <- raw_data@mzrange[2]
save(high_mz, file = "highest_mz.RData")

