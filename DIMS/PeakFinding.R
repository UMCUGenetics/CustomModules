#!/usr/bin/Rscript
# adapted from 4-peakFinding.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

filepath <- cmd_args[1]
breaks_filepath <- cmd_args[2] # location of breaks.fwhm.RData
resol <- as.numeric(cmd_args[3])
scripts_dir <- cmd_args[4]

thresh <- 2000

# if file extension is .txt, do nothing.
if (grepl(".txt$", filepath)) { stop("not run on txt file") }

# for debugging:
print(filepath)
print(breaks_filepath)
print(resol)
print(scripts_dir)
print(paste(scripts_dir, "AddOnFunctions/", sep=""))

# load in function scripts
# source(paste(scripts_dir, "AddOnFunctions/sourceDir.R", sep="/"))
# the sourceDir function no longer works
# sourceDir(paste(scripts_dir, "AddOnFunctions", sep="/"))
source(paste(scripts_dir, "AddOnFunctions/findPeaks.Gauss.HPC.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/searchMZRange.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/generateGaussian.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fitGaussian.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fitGaussianInit.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/getFwhm.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/getSD.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/optimizeGauss.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fit1Peak.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fit2peaks.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fit3peaks.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fit4peaks.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fitG.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fit2G.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fit3G.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/fit4G.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/getArea.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/getFitQuality.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/checkOverlap.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/sumCurves.R", sep=""))
source(paste(scripts_dir, "AddOnFunctions/isWithinXppm.R", sep=""))

load(breaks_filepath)
# load(filepath)

# for some reason, the repl_pattern_positive.RData and repl_pattern_negative.RData files are also read in.
if (!grepl("repl_pattern", filepath)) {
  # Load output of AverageTechReplicates for a sample
  sample_avgtechrepl <- get(load(filepath))
  if (grepl("_pos", filepath)) { scanmode = "positive" } else
    if (grepl("_neg", filepath)) { scanmode = "negative" }

  # for debugging:
  print(filepath)
  print(scanmode)

  # Initialize
  options(digits = 16)
  int.factor <- 1*10^5 # Number of x used to calc area under Gaussian (is not analytic)
  scale <- 2 # Initial value used to estimate scaling parameter
  width <- 1024
  height <- 768
  findPeaks.Gauss.HPC(sample_avgtechrepl, breaks.fwhm, int.factor, scale, resol, outdir, scanmode, FALSE, thresh, width, height)
}
