#!/usr/bin/Rscript
# Adapted from 1-generateBreaksFwhm.HPC.R

suppressPackageStartupMessages(library(mzR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

# Define parameters and variables
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")
filepath <- cmd_args[1]
outdir <- cmd_args[2]
trim <- as.numeric(cmd_args[3])  # 0.1
resol <- as.numeric(cmd_args[4])  # 140000
trimLeft=NULL
trimRight=NULL

# Open mzml file
Dat <- openMSfile(filepath)
hdr <- header(Dat)
pks <- spectra(Dat)


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

# Set trimLeft and trimRight (limits in terms of retention times)
trimLeft = round(hdr$retentionTime[2] + 0.5)
trimRight = round(hdr$retentionTime[length(hdr$retentionTime) - 1] - 0.5)

# Get the mass range m/z
lowMZ = round(min(hdr$lowMZ))
highMZ = round(max(hdr$highMZ))

# Create breaks/segmentation points within a range of m/z values with the aim
# to divide the m/z range into smaller sections/bins
nsegment = 2 * (highMZ - lowMZ)
segment = seq(from=lowMZ, to=highMZ, length.out = nsegment + 1)
breaks.fwhm = NULL
breaks.fwhm.avg = NULL

for (i in 1:nsegment) {
  startsegm <- segment[i]  # = current segment
  endsegm <- segment[i+1]
  resol.mz <- resol * (1/sqrt(2)^(log2(startsegm/200)))  # Calculates resolution at current segment
  fwhmsegm <- startsegm / resol.mz  # Calculates FWHM for peak in current segment
  breaks.fwhm <- c(breaks.fwhm, seq(from=(startsegm + fwhmsegm), to=endsegm, by=0.2 * fwhmsegm))  # Each segment is subdivided int smaller chunks (breaks) at intervals of 20% of the FWHM
  range = seq(from=(startsegm + fwhmsegm), to=endsegm, by=0.2 * fwhmsegm)  # Generates the same sequence of breakpoints as before but stores it in the variable range.
  deltaMZ = range[2]-range[1]
  breaks.fwhm.avg <- c(breaks.fwhm.avg, range + 0.5 * deltaMZ)  # This shifts the range values by half of deltaMZ to find the average m/z value within each bin
}

# Save various outputs in an R object
save(breaks.fwhm, breaks.fwhm.avg, trimLeft, trimRight, file = "breaks.fwhm.RData")
save(highMZ, file = "highest_mz.RData")
