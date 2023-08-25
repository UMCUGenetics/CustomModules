#!/usr/bin/Rscript
# adapted from 4-peakFinding.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

filepath        <- cmd_args[1]
breaks_filepath <- cmd_args[2] # location of breaks.fwhm.RData
resol           <- as.numeric(cmd_args[3])
scripts_dir     <- cmd_args[4]

thresh <- 2000

# load in function scripts
source(paste0(scripts_dir, "AddOnFunctions/findPeaks.Gauss.HPC.R"))
source(paste0(scripts_dir, "AddOnFunctions/searchMZRange.R"))
source(paste0(scripts_dir, "AddOnFunctions/generateGaussian.R"))
source(paste0(scripts_dir, "AddOnFunctions/fitGaussian.R"))
source(paste0(scripts_dir, "AddOnFunctions/fitGaussianInit.R"))
source(paste0(scripts_dir, "AddOnFunctions/getFwhm.R"))
source(paste0(scripts_dir, "AddOnFunctions/getSD.R"))
source(paste0(scripts_dir, "AddOnFunctions/optimizeGauss.R"))
source(paste0(scripts_dir, "AddOnFunctions/fit1Peak.R"))
source(paste0(scripts_dir, "AddOnFunctions/fit2peaks.R"))
source(paste0(scripts_dir, "AddOnFunctions/fit3peaks.R"))
source(paste0(scripts_dir, "AddOnFunctions/fit4peaks.R"))
source(paste0(scripts_dir, "AddOnFunctions/fitG.R"))
source(paste0(scripts_dir, "AddOnFunctions/fit2G.R"))
source(paste0(scripts_dir, "AddOnFunctions/fit3G.R"))
source(paste0(scripts_dir, "AddOnFunctions/fit4G.R"))
source(paste0(scripts_dir, "AddOnFunctions/getArea.R"))
source(paste0(scripts_dir, "AddOnFunctions/getFitQuality.R"))
source(paste0(scripts_dir, "AddOnFunctions/checkOverlap.R"))
source(paste0(scripts_dir, "AddOnFunctions/sumCurves.R"))
source(paste0(scripts_dir, "AddOnFunctions/isWithinXppm.R"))

load(breaks_filepath)

# Load output of AverageTechReplicates for a sample
sample_avgtechrepl <- get(load(filepath))
if (grepl("_pos", filepath)) { scanmode = "positive" } else
  if (grepl("_neg", filepath)) { scanmode = "negative" }

# Initialize
options(digits = 16)
int.factor <- 1*10^5 # Number of x used to calc area under Gaussian (is not analytic)
scale      <- 2 # Initial value used to estimate scaling parameter
width      <- 1024
height     <- 768

# run the findPeaks function
findPeaks.Gauss.HPC(sample_avgtechrepl, breaks.fwhm, int.factor, scale, resol, outdir, scanmode, FALSE, thresh, width, height)
