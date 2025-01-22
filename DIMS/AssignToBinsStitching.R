#!/usr/bin/Rscript
# Adapted from 2-DIMS.R

# Packages
suppressPackageStartupMessages({
  library(mzR)
  library(dplyr)
  library(stringr)
})

# Parameters and variables
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep = "")
filepath <- cmd_args[1]
breaks_filepath <- cmd_args[2]
resol <- as.numeric(cmd_args[3])
trim <- as.numeric(cmd_args[4])
dimsThresh <- 100
sampname <- sub('\\..*$', '', basename(filepath))
options(digits = 16)

# Open mzml and breaks files
Dat <- openMSfile(filepath)
hdr <- header(Dat)
pks <- spectra(Dat)
load(breaks_filepath)
bins <- rep(0, length(breaks.fwhm) - 1)
posBins <- bins
negBins <- bins
posRes <- NULL
negRes <- NULL

# FUNCTIONS: Adjust mz ranges
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


# FUNCTIONS: update filter string
Replace_FTMS <- function(hdr) {
  for (i in 1:length(hdr$filterString)) {
    hdr$filterString[i] <- str_replace(hdr$filterString[i], "(?<=\\[)\\d+", toString(hdr$scanWindowLowerLimit[i]))
  }
  return(hdr)
}

Replace_FTMShigh <- function(hdr) {
  for (i in 1:length(hdr$filterString)) {
    hdr$filterString[i] <- str_replace(hdr$filterString[i], "(?<=\\-)\\d+", toString(hdr$scanWindowUpperLimit[i]))
  }
  return(hdr)
}


# FUNCTIONS: Trim peaks by mz window
TrimPeaklistlower <- function(pks, hdr) {
  for (i in 1:length(pks)) {
    x <- which(pks[[i]][, 1] < hdr$scanWindowLowerLimit[i], arr.ind = TRUE)
    if (length(x) != 0) {pks[[i]] <- pks[[i]][-x, ]}
  }
  return(pks)
}

TrimPeaklistupper <- function(pks, hdr) {
  for (i in 1:length(pks)) {
    y <- which(pks[[i]][,1] < hdr$scanWindowUppperLimit[i], arr.ind = TRUE)
    if (length(y) != 0) {pks[[i]] <- pks[[i]][-y, ]}
  }
  return(pks)
}


# FUNCTION: Replaces peaksCount, lowMZ and HighMZ in hdr based on the new peak list pks
Replace_low_high_MZ_and_pkCount <- function(hdr, pks) {
  for (i in 1:length(pks)) {
    hdr$peaksCount[[i]] <- length(pks[[i]][, 1])
    hdr$lowMZ[[i]] <- min(pks[[i]][, 1])
    hdr$highMZ[[i]] <- max(pks[[i]][, 1])
  }
  return(hdr)
}


# Execute functions
hdr <- TrimlowerMZrange(hdr)
hdr <- TrimupperMZrange(hdr)
hdr <- Replace_FTMS(hdr)
hdr <- Replace_FTMShigh(hdr)
pks <- TrimPeaklistlower(pks, hdr)
pks <- TrimPeaklistupper(pks, hdr)
hdr <- Replace_low_high_MZ_and_pkCount(hdr, pks)

# Generate TIC table
TIC_intensity_persample <- cbind(round(hdr$retentionTime, 2), hdr$totIonCurrent)
colnames(TIC_intensity_persample) <- c("retention_time", "tic_intensity")
write.table(TIC_intensity_persample, file = paste0(sampname, "_TIC.txt"))

# Generate a matrix; y will contain a matrix with three columns: retention time,
# m/z, and intensity for all the scans combined
y <- do.call(rbind, lapply(1:length(pks), function(i) {
  rt <- hdr$retentionTime[i]
  scan_peaks <- pks[[i]]  
  cbind(rt, scan_peaks)  
 }))
colnames(y) <- c("time", "mz", "intensity")

# Get time values for positive and negative scans
posTimes <- hdr$retentionTime[hdr$polarity == 1]
negTimes <- hdr$retentionTime[hdr$polarity == 0]

# Filtering retention times
posTimes <- posTimes[posTimes > trimLeft & posTimes < trimRight]
negTimes <- negTimes[negTimes > trimLeft & negTimes < trimRight]

# Generating indices for scan modes
posInd <- which(y[, "time"] %in% posTimes)
negInd <- which(y[, "time"] %in% negTimes)  

# Separate scan modes into own matrices
posY <- y[posInd, ]
negY <- y[negInd, ]

# Divides the m/z values into discrete bins based on breaks.fwhm
yp <- cut(posY[, "mz"], breaks.fwhm, include.lowest = TRUE, right = TRUE, labels = FALSE) 
yn <- cut(negY[, "mz"], breaks.fwhm, include.lowest = TRUE, right = TRUE, labels = FALSE)

# Get the list of intensity values for each bin, and average the intensity values
# which are in the same bin
if (nrow(posY) > 0) {
  ap <- aggregate(posY[, "intensity"], list(yp), FUN = function(x) {
    if (is.na(mean(x[which(x > dimsThresh)]))) {  
      0 
    } else {
      mean(x[which(x > dimsThresh)])
    }
  })
  posBins[ap[, 1]] <- ap[, 2]
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

# Zero any values that are below the threshold
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
a <- breaks.fwhm.avg[-length(breaks.fwhm.avg)]  
b <- sprintf("%.5f", a)  
colnames(posRes) <- b
colnames(negRes) <- b

# Transpose matrices back to their original orientation, so that the m/z bins are in rows, and the samples are in columns
posResT <- t(posRes)
negResT <- t(negRes)

# Pack results in a list and save
peak_list <- list("pos" = posResT, "neg" = negResT, "breaksFwhm" = breaks.fwhm)
save(peak_list, file = paste0(sampname, ".RData"))
