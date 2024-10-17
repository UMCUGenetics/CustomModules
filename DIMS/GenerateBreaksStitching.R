#!/usr/bin/Rscript

## PACKAGES

suppressPackageStartupMessages(library(mzR))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))


## DEFINE PARAMETERS AND VARIABLES

cmd_args <- commandArgs(trailingOnly = TRUE)
#cmd_args <- c("/Users/nbuter/Documents/DIMS_data/DIMS_stitch_oefen_2/RES_20240827_UR_RUN7_002.mzML", "/Users/nbuter/Documents/DIMS_data/DIMS_stitch_oefen_2/", "0.1", "140000", "2", "1")  # NINA
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

filepath <- cmd_args[1]
outdir <- cmd_args[2]
trim <- as.numeric(cmd_args[3])  # 0.1
resol <- as.numeric(cmd_args[4])  # 140000

trimLeft=NULL
trimRight=NULL


## FUNCTIONS

# Trims the MZ range of scans in header file (lower scan window limit),  while
# leaving the scans with minimum values unchanged
TrimlowerMZrange <- function(hdr){
  for(i in 1:length(hdr$scanWindowLowerLimit)){
    #if (70%in%hdr$scanWindowLowerLimit[[i]]){print('Lowerlimit of current scan is 70')}
    if (hdr$scanWindowLowerLimit[[i]] == min(hdr$scanWindowLowerLimit)) {
      cat('Lowerlimit has reached:', hdr$scanWindowLowerLimit[[i]])
      }
    else {hdr$scanWindowLowerLimit[[i]] <- hdr$scanWindowLowerLimit[[i]] + 5}
  }
  return(hdr)
}


# Trims the MZ range of scans in header file (upper scan window limit), while
# leaving the scans with maximum values unchanged
TrimupperMZrange <- function(hdr){
  for(i in 1:length(hdr$scanWindowUpperLimit)){
    #if (1280%in%hdr$scanWindowUpperLimit[[i]]){print('Upperlimit of current scan is 1280')}
    if (hdr$scanWindowUpperLimit[[i]] == max(hdr$scanWindowUpperLimit)) {
      cat('Upperlimit has reached:', hdr$scanWindowUpperLimit[[i]])
      }
    else {hdr$scanWindowUpperLimit[[i]] <- hdr$scanWindowUpperLimit[[i]] - 5}
  }
  return(hdr)
}


# Replace lower mz in hdr filterString column with correct mz created with
# TrimlowerMZrange
Replace_FTMS <- function(hdr){
  for(i in 1:length(hdr$filterString)){
    hdr$filterString[i] <- str_replace(hdr$filterString[i],"(?<=\\[)\\d+", toString(hdr$scanWindowLowerLimit[i]))
  }
  return(hdr)
}


# Replace upper mz in hdr filterString column with correct mz created with
# TrimupperMZrange
Replace_FTMShigh <- function(hdr){
  for(i in 1:length(hdr$filterString)) {
    hdr$filterString[i] <- str_replace(hdr$filterString[i],"(?<=\\-)\\d+", toString(hdr$scanWindowUpperLimit[i]))
  }
  return(hdr)
}


# Removes all m/z peaks in the peak list pks that are below the lower limit of
# the m/z window for each scan
TrimPeaklistlower <- function(pks, hdr){
  for (i in 1:length(pks)) {
    x <- which(pks[[i]][,1] < hdr$scanWindowLowerLimit[i], arr.ind = TRUE)
    if (length(x)!=0) {pks[[i]] <- pks[[i]][-x,]}
  }
  return(pks)
}


# Removes all m/z peaks in the peak list pks that are above the upper limit of
# the m/z window for each scan
TrimPeaklistupper <- function(pks, hdr){
  for (i in 1:length(pks)) {
    y <- which(pks[[i]][,1] < hdr$scanWindowUppperLimit[i], arr.ind = TRUE)
    if (length(y)!=0) {pks[[i]] <- pks[[i]][-y,]}
  }
  return(pks)
}


# Replaces peaksCount, lowMZ and HighMZ in hdr based on the new peak list pks
Replace_low_high_MZ_and_pkCount <- function(hdr, pks){
  for (i in 1:length(pks)) {
    hdr$peaksCount[[i]] <- length(pks[[i]][,1])
    hdr$lowMZ[[i]] <- min(pks[[i]][,1])
    hdr$highMZ[[i]] <- max(pks[[i]][,1])
  }
  return(hdr)
}


## OPEN mzML FILE AND EXECUTE FUNCTIONS (NEEDED FOR OBTAINING MASS RANGE FOR CREATING BREAKS)


Dat <- openMSfile(filepath)
hdr <- header(Dat)
pks <- spectra(Dat)
hdr <- TrimlowerMZrange(hdr)
hdr <- TrimupperMZrange(hdr)
hdr <- Replace_FTMS(hdr)
hdr <- Replace_FTMShigh(hdr)
pks <- TrimPeaklistlower(pks, hdr)
pks <- TrimPeaklistupper(pks, hdr)
hdr <- Replace_low_high_MZ_and_pkCount(hdr, pks)


## GET MASS RANGE, CREATE, AND SAVE THE BREAKS FILE

# Set trimLeft and trimRight (limits in terms of retention times)
trimLeft = round(hdr$retentionTime[2] + 0.5)
trimRight = round(hdr$retentionTime[length(hdr$retentionTime) - 1] - 0.5)
cat(paste("\ntrimLeft", trimLeft, sep=" "))
cat(paste("\ntrimRight", trimRight, sep=" "))


# Get the mass range m/z
lowMZ = round(min(hdr$lowMZ))
highMZ = round(max(hdr$highMZ))


# Create breaks/segmentation points within a range of m/z values with the aim
# to divide the m/z range into smaller sections/bins
nsegment = 2 * (highMZ - lowMZ)
segment = seq(from=lowMZ, to=highMZ, length.out = nsegment + 1)  # Creates sequence of nsegments + 1 equally spaced points
breaks.fwhm = NULL
breaks.fwhm.avg = NULL


for (i in 1:nsegment) {
  startsegm <- segment[i]  # = current segment
  endsegm <- segment[i+1]
  resol.mz <- resol * (1/sqrt(2)^(log2(startsegm/200)))  # Calculates resolution at current segment
  fwhmsegm <- startsegm / resol.mz  # Calculates FWHM for peak in current segment
  breaks.fwhm <- c(breaks.fwhm, seq(from=(startsegm + fwhmsegm), to=endsegm, by=0.2 * fwhmsegm))  # Each segment is subdivided int smaller chunks (breaks) at intervals of 20% of the FWHM
  #breaks.fwhm <- c(breaks.fwhm, seq(from=(startsegm), to=endsegm, by=0.2*fwhmsegm))
  # average the m/z instead of start value (?)
  range = seq(from=(startsegm + fwhmsegm), to=endsegm, by=0.2 * fwhmsegm)  # Generates the same sequence of breakpoints as before but stores it in the variable range.
  deltaMZ = range[2]-range[1]
  breaks.fwhm.avg <- c(breaks.fwhm.avg, range + 0.5 * deltaMZ)  # This shifts the range values by half of deltaMZ to find the average m/z value within each bin
}


# Save various outputs in an R object
save(breaks.fwhm, breaks.fwhm.avg, trimLeft, trimRight, file = "breaks.fwhm.RData")
save(highMZ, file = "highest_mz.RData")