#!/usr/bin/Rscript
## adapted from 2-DIMS.R

# load required packages
suppressPackageStartupMessages(library("xcms"))

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

filepath        <- cmd_args[1]
breaks_filepath <- cmd_args[2]
resol           <- as.numeric(cmd_args[3])
trim            <- 0.1
dims_thresh     <- 100

# get sample name
sample_name <- sub("\\..*$", "", basename(filepath))

options(digits = 16)

# Initialize
int_factor <- 1 * 10^5 # Number used to calculate area under Gaussian curve
scale      <- 2 # Initial value used to estimate scaling parameter
width      <- 1024
height     <- 768
trim_left       <- NULL
trim_right      <- NULL
breaks_fwhm     <- NULL
breaks_fwhm_avg <- NULL
bins            <- NULL
pos_results     <- NULL
neg_results     <- NULL

# read in the data for 1 sample
raw_data <- suppressMessages(xcmsRaw(filepath))

# for TIC plots: prepare txt files with data for plots
tic_intensity_persample <- cbind(round(raw_data@scantime, 2), raw_data@tic)
colnames(tic_intensity_persample) <- c("retention_time", "tic_intensity")
write.table(tic_intensity_persample, file = paste0("./", sample_name, "_TIC.txt"))

# load breaks_fwhm
load(breaks_filepath)

# Create empty placeholders for later use
bins <- rep(0, length(breaks_fwhm) - 1)
pos_bins <- bins
neg_bins <- bins

# Generate a matrix
raw_data_matrix <- rawMat(raw_data)

# Get time values for positive and negative scans
pos_times <- raw_data@scantime[raw_data@polarity == "positive"]
neg_times <- raw_data@scantime[raw_data@polarity == "negative"]
# Select scans between trim_left and trim_right
pos_times <- pos_times[pos_times > trim_left & pos_times < trim_right]
neg_times <- neg_times[neg_times > trim_left & neg_times < trim_right]

# Generate an index with which to select values for each mode
pos_index <- which(raw_data_matrix[, "time"] %in% pos_times)
neg_index <- which(raw_data_matrix[, "time"] %in% neg_times)
# Separate each mode into its own matrix
pos_raw_data_matrix <- raw_data_matrix[pos_index, ]
neg_raw_data_matrix <- raw_data_matrix[neg_index, ]

# Get index for binning intensity values
bin_indices_pos <- cut(pos_raw_data_matrix[, "mz"], breaks_fwhm,
                       include.lowest = TRUE, right = TRUE, labels = FALSE)
bin_indices_neg <- cut(neg_raw_data_matrix[, "mz"], breaks_fwhm,
                       include.lowest = TRUE, right = TRUE, labels = FALSE)

# Get the list of intensity values for each bin, and add the
# intensity values which are in the same bin
if (nrow(pos_raw_data_matrix) > 0) {
  # set NA in intensities to zero
  pos_raw_data_matrix[is.na(pos_raw_data_matrix[, "intensity"]), "intensity"] <- 0
  # use only values above dims_thresh
  pos_intensity_above_threshold <- pos_raw_data_matrix[which(pos_raw_data_matrix[, "intensity"] > dims_thresh), "intensity"]
  # aggregate intensities, calculate mean
  aggr_int_pos <- stats::aggregate(pos_intensity_above_threshold, list(bin_indices_pos), mean)
  pos_bins[aggr_int_pos[, 1]] <- aggr_int_pos[, 2]
}
if (nrow(neg_raw_data_matrix) > 0) {
  # set NA in intensities to zero
  neg_raw_data_matrix[is.na(neg_raw_data_matrix[, "intensity"]), "intensity"] <- 0
  # use only values above dims_thresh
  neg_intensity_above_threshold <- neg_raw_data_matrix[which(neg_raw_data_matrix[, "intensity"] > dims_thresh), "intensity"]
  # aggregate intensities, calculate mean
  aggr_int_neg <- stats::aggregate(neg_intensity_above_threshold, list(bin_indices_neg), mean)
  neg_bins[aggr_int_neg[, 1]] <- aggr_int_neg[, 2]
}

# Zero any values that are below the threshold
pos_bins[pos_bins < dims_thresh] <- 0
neg_bins[neg_bins < dims_thresh] <- 0

pos_results <- cbind(pos_results, pos_bins)
neg_results <- cbind(neg_results, neg_bins)

# transpose
pos_results_transpose <- t(pos_results)
neg_results_transpose <- t(neg_results)

# Add file names as row names
rownames(pos_results_transpose) <- sample_name
rownames(neg_results_transpose) <- sample_name

# delete the last value of breaks_fwhm_avg to match dimensions of pos_results and neg_results
breaks_fwhm_avg_minus1 <- breaks_fwhm_avg[-length(breaks_fwhm_avg)]
# Format as string and show precision of float to 5 digits
breaks_fwhm_avg_minus1 <- sprintf("%.5f", breaks_fwhm_avg_minus1)

# Use this as the column names
colnames(pos_results_transpose) <- breaks_fwhm_avg_minus1
colnames(neg_results_transpose) <- breaks_fwhm_avg_minus1

# transpose back
pos_results_final <- t(pos_results_transpose)
neg_results_final <- t(neg_results_transpose)

peak_list <- list("pos" = pos_results_final, "neg" = neg_results_final, "breaksFwhm" = breaks_fwhm)

save(peak_list, file = paste0("./", sample_name, ".RData"))