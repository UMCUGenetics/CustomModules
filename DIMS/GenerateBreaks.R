#### GenerateBreaks.R ####
## adapted from 1-generateBreaksFwhm.HPC.R ##
#!/usr/bin/Rscript

# .libPaths(new = "/hpc/local/CentOS7/dbg_mz/R_libs/3.2.2")

# load required package 
suppressPackageStartupMessages(library("xcms"))

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

filepath <- cmd_args[1] # 1 of the mzML files
outdir   <- cmd_args[2] 
trim     <- as.numeric(cmd_args[3]) # 0.1
resol    <- as.numeric(cmd_args[4]) # 140000

# initialize
trimLeft        = NULL
trimRight       = NULL
breaks.fwhm     = NULL
breaks.fwhm.avg = NULL
bins            = NULL
posRes          = NULL
negRes          = NULL

# read in mzML file
raw_data <- suppressMessages(xcmsRaw(filepath))

# trim scans at the start and end
trimLeft  = round(raw_data@scantime[length(raw_data@scantime)*trim])
trimRight = round(raw_data@scantime[length(raw_data@scantime)*(1-trim)])

# Mass range m/z
lowMZ = raw_data@mzrange[1]
highMZ = raw_data@mzrange[2]

# determine number of segments (bins)
nsegment = 2*(highMZ-lowMZ)
segment  = seq(from=lowMZ, to=highMZ, length.out=nsegment+1)

# determine start and end of each bin.
for (i in 1:nsegment) {
  startsegm   <- segment[i]
  endsegm     <- segment[i+1]
  resol.mz    <- resol*(1/sqrt(2)^(log2(startsegm/200)))
  fwhmsegm    <- startsegm/resol.mz
  breaks.fwhm <- c(breaks.fwhm, seq(from=(startsegm + fwhmsegm), to=endsegm, by=0.2*fwhmsegm))
  # average the m/z instead of start value
  range = seq(from=(startsegm + fwhmsegm), to=endsegm, by=0.2*fwhmsegm)
  deltaMZ = range[2] - range[1]
  breaks.fwhm.avg <- c(breaks.fwhm.avg, range + 0.5*deltaMZ)
}

save(breaks.fwhm, breaks.fwhm.avg, trimLeft, trimRight, file=paste(outdir, "breaks.fwhm.RData", sep="/"))

