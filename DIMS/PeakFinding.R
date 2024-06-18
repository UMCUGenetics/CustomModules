## adapted from 4-peakFinding.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

sample_file <- cmd_args[1]
breaks_file <- cmd_args[2]
resol <- as.numeric(cmd_args[3])
scripts_dir <- cmd_args[4]
thresh <- 2000
outdir <- "./"

# load in function scripts
source(paste0(scripts_dir, "do_peakfinding.R"))
source(paste0(scripts_dir, "check_overlap.R"))
source(paste0(scripts_dir, "search_mzrange.R"))
source(paste0(scripts_dir, "fit_optim.R"))
source(paste0(scripts_dir, "fit_gaussian.R"))
source(paste0(scripts_dir, "fit_init.R"))
source(paste0(scripts_dir, "get_fwhm.R"))
source(paste0(scripts_dir, "get_stdev.R"))
source(paste0(scripts_dir, "optimize_gaussfit.R"))
source(paste0(scripts_dir, "fit_peaks.R"))
source(paste0(scripts_dir, "fit_gaussians.R"))
source(paste0(scripts_dir, "estimate_area.R"))
source(paste0(scripts_dir, "get_fit_quality.R"))
source(paste0(scripts_dir, "check_overlap.R"))
source(paste0(scripts_dir, "sum_curves.R"))
source(paste0(scripts_dir, "within_ppm.R"))

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
# Number used to calculate area under Gaussian curve
int_factor <- 1 * 10^5 
# Initial value used to estimate scaling parameter
scale <- 2
width <- 1024
height <- 768

# run the findPeaks function

# do_peakfinding(sample_avgtechrepl, breaks_fwhm, int_factor, scale, resol, outdir, scanmode, FALSE, thresh, width, height)
do_peakfinding(sample_avgtechrepl, int_factor, scale, resol, outdir, scanmode, FALSE, thresh, width, height)
