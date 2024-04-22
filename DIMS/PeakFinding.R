#!/usr/bin/Rscript
# adapted from 4-peakFinding.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

sample_file <- cmd_args[1]
breaks_file <- cmd_args[2]
resol <- as.numeric(cmd_args[3])
scripts_dir <- cmd_args[4]
thresh <- 2000
outdir <- "./"

# load in function scripts
source(paste0(scripts_dir, "findPeaks.Gauss.HPC.R"))
source(paste0(scripts_dir, "searchMZRange.R"))
source(paste0(scripts_dir, "generateGaussian.R"))
source(paste0(scripts_dir, "fitGaussian.R"))
source(paste0(scripts_dir, "fitGaussianInit.R"))
source(paste0(scripts_dir, "getFwhm.R"))
source(paste0(scripts_dir, "getSD.R"))
source(paste0(scripts_dir, "optimizeGauss.R"))
source(paste0(scripts_dir, "fit1Peak.R"))
source(paste0(scripts_dir, "fit2peaks.R"))
source(paste0(scripts_dir, "fit3peaks.R"))
source(paste0(scripts_dir, "fit4peaks.R"))
source(paste0(scripts_dir, "fitG.R"))
source(paste0(scripts_dir, "fit2G.R"))
source(paste0(scripts_dir, "fit3G.R"))
source(paste0(scripts_dir, "fit4G.R"))
source(paste0(scripts_dir, "getArea.R"))
source(paste0(scripts_dir, "getFitQuality.R"))
source(paste0(scripts_dir, "checkOverlap.R"))
source(paste0(scripts_dir, "sumCurves.R"))
source(paste0(scripts_dir, "isWithinXppm.R"))

load(breaks_file)

# Load output of AverageTechReplicates for a sample
sample_avgtechrepl <- get(load(sample_file))
if (grepl("_pos", sample_file)) {
  scanmode <- "positive"
} else if (grepl("_neg", sample_file)) {
  scanmode <- "negative"
}

# Initialize
options(digits = 16)
int_factor <- 1 * 10^5 # Number used to calculate area under Gaussian curve
scale <- 2 # Initial value used to estimate scaling parameter
width <- 1024
height <- 768

# run the findPeaks function
print(head(sample_avgtechrepl))
print(head(breaks_fwhm))
print(int_factor)
print(scale)
print(resol)
print(outdir)
print(scanmode)
print(thresh)
print(width)
print(height)

do_peakfinding(sample_avgtechrepl, breaks_fwhm, int_factor, scale, resol, outdir, scanmode, FALSE, thresh, width, height)
