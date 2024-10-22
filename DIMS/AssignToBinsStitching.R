#!/usr/bin/Rscript


## PACKAGES

suppressPackageStartupMessages(library("mzR"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))


# DEFINE PARAMETERS AND VARIABLES

cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep = "")

filepath <- cmd_args[1]
breaks_filepath <- cmd_args[2]
resol <- as.numeric(cmd_args[3])
trim <- as.numeric(cmd_args[4])
dimsThresh <- 100


# Print sample names and for use later on
sampname <- sub('\\..*$', '', basename(filepath))

options(digits = 16)

# Define variables
trimLeft <- NULL
trimRight <- NULL
breaks.fwhm <- NULL
breaks.fwhm.avg <- NULL
bins <- NULL
posRes <- NULL
negRes <- NULL


## FUNCTIONS FOR STITCHING

# Trims the MZ range of scans in header file (lower scan window limit),  while
# leaving the scans with minimum values unchanged
TrimlowerMZrange <- function(hdr) {
  cat("\n trim lower MZ range:")
  for (i in 1:length(hdr$scanWindowLowerLimit)) {
    if (hdr$scanWindowLowerLimit[[i]] == min(hdr$scanWindowLowerLimit)) {
      cat(',', hdr$scanWindowLowerLimit[[i]])
      }
    else {hdr$scanWindowLowerLimit[[i]] <- hdr$scanWindowLowerLimit[[i]] + 5}
  }
  return(hdr)
}


# Trims the MZ range of scans in header file (upper scan window limit), while
# leaving the scans with maximum values unchanged
TrimupperMZrange <- function(hdr) {
  cat("\n trim upper MZ range:")
  for (i in 1:length(hdr$scanWindowUpperLimit)) {
    if (hdr$scanWindowUpperLimit[[i]] == max(hdr$scanWindowUpperLimit)) {
      cat(',', hdr$scanWindowUpperLimit[[i]])
      }
    else {hdr$scanWindowUpperLimit[[i]] <- hdr$scanWindowUpperLimit[[i]] - 5}
  }
  return(hdr)
}


# Replace lower mz in hdr filterString column with correct mz created with
# TrimlowerMZrange
Replace_FTMS <- function(hdr) {
  for (i in 1:length(hdr$filterString)) {
    hdr$filterString[i] <- str_replace(hdr$filterString[i], "(?<=\\[)\\d+", toString(hdr$scanWindowLowerLimit[i]))
  }
  return(hdr)
}


# Replace upper mz in hdr filterString column with correct mz created with
# TrimupperMZrange
Replace_FTMShigh <- function(hdr) {
  for (i in 1:length(hdr$filterString)) {
    hdr$filterString[i] <- str_replace(hdr$filterString[i], "(?<=\\-)\\d+", toString(hdr$scanWindowUpperLimit[i]))
  }
  return(hdr)
}


# Removes all m/z peaks in the peak list pks that are below the lower limit of
# the m/z window for each scan
TrimPeaklistlower <- function(pks, hdr) {
  for (i in 1:length(pks)) {
    x <- which(pks[[i]][, 1] < hdr$scanWindowLowerLimit[i], arr.ind = TRUE)
    if (length(x) != 0) {pks[[i]] <- pks[[i]][-x, ]}
  }
  return(pks)
}


# Removes all m/z peaks in the peak list pks that are above the upper limit of
# the m/z window for each scan
TrimPeaklistupper <- function(pks, hdr) {
  for (i in 1:length(pks)) {
    y <- which(pks[[i]][,1] < hdr$scanWindowUppperLimit[i], arr.ind = TRUE)
    if (length(y) != 0) {pks[[i]] <- pks[[i]][-y, ]}
  }
  return(pks)
}


# Replaces peaksCount, lowMZ and HighMZ in hdr based on the new peak list pks
Replace_low_high_MZ_and_pkCount <- function(hdr, pks) {
  for (i in 1:length(pks)) {
    hdr$peaksCount[[i]] <- length(pks[[i]][, 1])
    hdr$lowMZ[[i]] <- min(pks[[i]][, 1])
    hdr$highMZ[[i]] <- max(pks[[i]][, 1])
  }
  return(hdr)
}


## OPEN mzML FILE AND EXECUTE FUNCTIONS IF STITCH DATA

# Open files
Dat <- openMSfile(filepath)
hdr <- header(Dat)
pks <- spectra(Dat)

# Execute functions
hdr <- TrimlowerMZrange(hdr)
hdr <- TrimupperMZrange(hdr)
hdr <- Replace_FTMS(hdr)
hdr <- Replace_FTMShigh(hdr)
pks <- TrimPeaklistlower(pks, hdr)
pks <- TrimPeaklistupper(pks, hdr)
hdr <- Replace_low_high_MZ_and_pkCount(hdr, pks)

# Creates table for TIC plots
TIC_intensity_persample <- cbind(round(hdr$retentionTime, 2), hdr$totIonCurrent)
colnames(TIC_intensity_persample) <- c("retentionTime", "TIC")
write.table(TIC_intensity_persample, file = paste0(sampname, "_TIC.txt"))


## ASSIGN THE POSITIVE AND NEGATIVE ION MODES' PEAKS (INTENSITIES) TO BINS

# Load in breaks.fwhm.RData
load(breaks_filepath)

# Create empty placeholders for later use
bins <- rep(0, length(breaks.fwhm) - 1)

# Generate a matrix; y will contain a matrix with three columns: retention time,
# m/z, and intensity for all the scans combined
times <- hdr$retentionTime
y <- NA
t <- NULL
for (i in 1:length(pks)) {
  t <- times[i]
  this_scan <- pks[[i]]
  l <- length(this_scan[, 1])  # Stores the number of peaks (rows) in this scan
  time <- rep(t, l)  # Creates a vector time of length l, where each element is the retention time t
  this_scan <- cbind(time, this_scan) # Combines the retention time time with the peak data (this_scan), so now each row contains: retention time, m/z, and intensity
  colnames(this_scan) <- c("time", "mz", "intensity")
  if (i == 1) {y <- this_scan} else {y <- rbind(y, this_scan)}
}


# This might be a more efficient method to combine the matrices from pks with
# retention times from hdr than the for-loop above (NINA)
# y2 <- do.call(rbind, lapply(1:length(pks), function(i) {
#   rt <- times[i] # Get retention time from current scan
#   scan_peaks <- pks[[i]]  # Get peak data from current scan
#   cbind(rt, scan_peaks)  # Combine
# }))


# Set the column names (NINA)
# colnames(y2) <- c("time", "mz", "intensity")


# Get time values for positive and negative scans
posTimes <- hdr$retentionTime[hdr$polarity == 1]
negTimes <- hdr$retentionTime[hdr$polarity == 0]


# Set trimLeft and trimRight (limits in terms of retention times) stitching vs.
# non-stitching data (NOT NECESSARY - NINA)
#trimLeft <- round(hdr$retentionTime[2] + 0.5)
#trimRight <- round(hdr$retentionTime[length(hdr$retentionTime)-1] - 0.5)


# Filtering retention times
posTimes <- posTimes[posTimes > trimLeft & posTimes < trimRight]
negTimes <- negTimes[negTimes > trimLeft & negTimes < trimRight]


# Generating index for selections (positive and negative ion mode)
posInd <- which(y[, "time"] %in% posTimes)  # Returns the indices of rows where retention times in y ("time") match those found in posTimes
negInd <- which(y[, "time"] %in% negTimes)  # Same but for the negative mode (negTimes)


# Separate scans of each mode (positive and negative ion modes) into own matrix
posY <- y[posInd, ]
negY <- y[negInd, ]


# Divides the m/z values into discrete bins based on breaks.fwhm
## This doesn't round the value for mz - is this an issue?
yp <- cut(posY[, "mz"], breaks.fwhm, include.lowest = TRUE, right = TRUE, labels = FALSE) # This will contain the bin index for each m/z value in posY, indicating which segment it belongs to based on the FWHM bins
yn <- cut(negY[, "mz"], breaks.fwhm, include.lowest = TRUE, right = TRUE, labels = FALSE)


# Empty the bins
posBins <- bins
negBins <- bins


# Get the list of intensity values for each bin, and average the intensity values
# which are in the same bin
if (nrow(posY) > 0) {
  ap <- aggregate(posY[, "intensity"], list(yp), FUN = function(x) {
    if (is.na(mean(x[which(x > dimsThresh)]))) {  # Filters intensities by keeping only values greater than dimsThresh
      0  # Function returns zero when all intensities are below threshold or bin has no values left after filtering
    } else {
      mean(x[which(x > dimsThresh)])  # Calculates the mean intensity of the remaining values for each m/z bin
    }
  })
  posBins[ap[, 1]] <- ap[, 2] # Replace existing values with new mean intensities

}


if (nrow(negY) > 0) {
  an <- aggregate(negY[, "intensity"], list(yn), FUN = function(x) {
    if (is.na(mean(x[which(x > dimsThresh)]))) {
      0
    } else {
      mean(x[which(x > dimsThresh)])
    }
  })
  negBins[an[, 1]] <- an[, 2]
}


# Zero any values that are below the threshold (low-intensity values)
posBins[posBins < dimsThresh] <- 0
negBins[negBins < dimsThresh] <- 0


# Creates results matrices
posRes <- cbind(posRes, posBins)
negRes <- cbind(negRes, negBins)


# Transpose the result matrices
posRes <- t(posRes)
negRes <- t(negRes)


# Add sample names as row names (one row per sample)
rownames(posRes) <- sampname
rownames(negRes) <- sampname


# Adjust breaks.fwhm.avg and create column names
a <- breaks.fwhm.avg[-length(breaks.fwhm.avg)]  # Removes the last element from breaks.fwhm.avg # + 0.5*deltaMZ (?)
b <- sprintf("%.5f", a)  # Formats the remaining values as strings with 5 decimal places of precision
colnames(posRes) <- b  # The result is used as column names, where each column corresponds to an m/z bin with a specific boundary
colnames(negRes) <- b


# Transpose matrices back to their original orientation, so that the m/z bins are in rows, and the samples are in columns
posResT <- t(posRes)
negResT <- t(negRes)


# Omit rows with only zeros (was commented out?)
#  sumsp <- apply(posResT,1,sum)
#  posResT.nonzero <- posResT[(sumsp != 0), ] # <=============================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  sums <- apply(negResT,1,sum)
#  negResT.nonzero <- negResT[(sums != 0), ] # <=============================================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# Pack results in a list and save
peak_list <- list("pos" = posResT, "neg" = negResT, "breaksFwhm" = breaks.fwhm)
save(peak_list, file = paste0(sampname, ".RData"))
