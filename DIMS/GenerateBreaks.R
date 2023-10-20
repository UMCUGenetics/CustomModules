#### GenerateBreaks.R ####
## adapted from 1-generateBreaksFwhm.HPC.R ##
#!/usr/bin/Rscript

# load required package 
suppressPackageStartupMessages(library("xcms"))

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)

filepath <- cmd_args[1] # 1 of the mzML files
outdir   <- cmd_args[2] 
trim     <- as.numeric(cmd_args[3]) # 0.1
resol    <- as.numeric(cmd_args[4]) # 140000

# initialize
trimLeft        <- NULL
trimRight       <- NULL
breaks_fwhm     <- NULL
breaks_fwhm_avg <- NULL
bins            <- NULL
posRes          <- NULL
negRes          <- NULL

# read in mzML file
raw_data <- suppressMessages(xcmsRaw(filepath))

# trim (remove) scans at the start and end
trimLeft  <- round(raw_data@scantime[length(raw_data@scantime)*trim])
trimRight <- round(raw_data@scantime[length(raw_data@scantime)*(1-trim)])

# Mass range m/z
lowMZ  <- raw_data@mzrange[1]
highMZ <- raw_data@mzrange[2]

# determine number of segments (bins)
nsegment <- 2*(highMZ-lowMZ)
segment  <- seq(from=lowMZ, to=highMZ, length.out=nsegment+1)

# determine start and end of each bin.
for (i in 1:nsegment) {
  startsegm   <- segment[i]
  endsegm     <- segment[i+1]
  resol_mz    <- resol*(1/sqrt(2)^(log2(startsegm/200)))
  fwhmsegm    <- startsegm/resol_mz
  breaks_fwhm <- c(breaks_fwhm, seq(from=(startsegm + fwhmsegm), to=endsegm, by=0.2*fwhmsegm))
  # average the m/z instead of start value
  range <- seq(from=(startsegm + fwhmsegm), to=endsegm, by=0.2*fwhmsegm)
  deltaMZ <- range[2] - range[1]
  breaks_fwhm_avg <- c(breaks_fwhm_avg, range + 0.5*deltaMZ)
}

# generate output file
save(breaks_fwhm, breaks_fwhm_avg, trimLeft, trimRight, file="./breaks.fwhm.RData")

