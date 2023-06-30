#!/usr/bin/Rscript

#.libPaths(new = "/hpc/local/CentOS7/dbg_mz/R_libs/3.6.2")

# load required packages 
suppressPackageStartupMessages(library("xcms"))

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

filepath <- cmd_args[1] # location of mzML file
breaks_filepath <- cmd_args[2] # location of breaks.fwhm.RData
resol <- as.numeric(cmd_args[3]) # 140000
trim <- 0.1
dimsThresh <- 100
print(filepath)
print(breaks_filepath)
print(resol)

#sampname <- sub('\\..*$', '', basename(filepath))
#cat(paste0("\n", sampname))
sampname <- "AssignToBins"

#suppressPackageStartupMessages(library("Cairo"))
options(digits=16)

### process one sample at a time and find peaks FOR BOTH SCAN MODES! #
int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic)
scale=2 # Initial value used to estimate scaling parameter
width=1024
height=768

# Aggregate
trimLeft=NULL
trimRight=NULL
breaks.fwhm=NULL
breaks.fwhm.avg=NULL
bins=NULL
pos_results=NULL
neg_results=NULL

# read in the data for 1 sample
raw_data <- suppressMessages(xcmsRaw(filepath))

# load breaks.fwhm
load(breaks_filepath)

# Create empty placeholders for later use
bins <- rep(0, length(breaks.fwhm) - 1)
pos_bins <- bins
neg_bins <- bins

# Generate a matrix
raw_data_matrix <- rawMat(raw_data)

# Get time values for positive and negative scans
pos_times <- raw_data@scantime[raw_data@polarity == "positive"]
neg_times <- raw_data@scantime[raw_data@polarity == "negative"]
# Select scans between trimLeft and trimRight
pos_times <- pos_times[pos_times > trimLeft & pos_times < trimRight]
neg_times <- neg_times[neg_times > trimLeft & neg_times < trimRight]

# Generate an index with which to select values for each mode
pos_index <- which(raw_data_matrix[ ,"time"] %in% pos_times)
neg_index <- which(raw_data_matrix[ ,"time"] %in% neg_times)
# Separate each mode into its own matrix
pos_raw_data_matrix <- raw_data_matrix[pos_index, ]
neg_raw_data_matrix <- raw_data_matrix[neg_index, ]

# Get index for binning intensity values
bin_indices_pos <- cut(pos_raw_data_matrix[ ,"mz"], breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)
bin_indices_neg <- cut(neg_raw_data_matrix[ ,"mz"], breaks.fwhm, include.lowest=TRUE, right=TRUE, labels=FALSE)

# Get the list of intensity values for each bin, and add the
# intensity values which are in the same bin
if (nrow(pos_raw_data_matrix) > 0) {
  # set NA in intensities to zero
  pos_raw_data_matrix[is.na(pos_raw_data_matrix[,"intensity"]), "intensity"] <- 0
  # use only values above dimsThresh
  pos_intensity_above_threshold <- pos_raw_data_matrix[which(pos_raw_data_matrix[ ,"intensity"] > dimsThresh), "intensity"]
  # aggregate intensities, calculate mean
  aggr_int_pos <- stats::aggregate(pos_intensity_above_threshold, list(bin_indices_pos), mean)
  pos_bins[aggr_int_pos[ ,1]] <- aggr_int_pos[ ,2]
}
if (nrow(neg_raw_data_matrix) > 0) {
  # set NA in intensities to zero
  neg_raw_data_matrix[is.na(neg_raw_data_matrix[,"intensity"]), "intensity"] <- 0
  # use only values above dimsThresh
  neg_intensity_above_threshold <- neg_raw_data_matrix[which(neg_raw_data_matrix[ ,"intensity"] > dimsThresh), "intensity"]
  # aggregate intensities, calculate mean
  aggr_int_neg <- stats::aggregate(neg_intensity_above_threshold, list(bin_indices_neg), mean)
  neg_bins[aggr_int_neg[ ,1]] <- aggr_int_neg[ ,2]
}

# Zero any values that are below the threshold
pos_bins[pos_bins < dimsThresh] <- 0
neg_bins[neg_bins < dimsThresh] <- 0

pos_results = cbind(pos_results, pos_bins)
neg_results = cbind(neg_results, neg_bins)

# transpose
pos_results_transpose = t(pos_results)
neg_results_transpose = t(neg_results)

# Add file names as row names
rownames(pos_results_transpose) = sampname
rownames(neg_results_transpose) = sampname

# delete the last value of breaks.fwhm.avg to match dimensions of pos_results and neg_results
breaks.fwhm.avg.minus1 <- breaks.fwhm.avg[-length(breaks.fwhm.avg)]
# Format as string and show precision of float to 5 digits
breaks.fwhm.avg.minus1 <- sprintf("%.5f", breaks.fwhm.avg.minus1)

# Use this as the column names
colnames(pos_results_transpose) <- breaks.fwhm.avg.minus1
colnames(neg_results_transpose) <- breaks.fwhm.avg.minus1

# transpose back
pos_results_final <- t(pos_results_transpose)
neg_results_final <- t(neg_results_transpose)

pklist <- list("pos"=pos_results_final, "neg"=neg_results_final, "breaksFwhm"=breaks.fwhm)

save(pklist, file=paste("./", sampname, ".RData", sep=""))
