#!/usr/bin/Rscript
# adapted from 9-runFillMissing.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

peakgrouplist_file <- cmd_args[1]
scripts_dir <- cmd_args[2]
thresh <- as.numeric(cmd_args[3])
resol <- as.numeric(cmd_args[4])
ppm <- as.numeric(cmd_args[5])
outdir <- "./"

# load in function scripts
source(paste0(scripts_dir, "replaceZeros.R"))
source(paste0(scripts_dir, "generateGaussian.R"))
source(paste0(scripts_dir, "getFwhm.R"))
source(paste0(scripts_dir, "getSD.R"))
source(paste0(scripts_dir, "getArea.R"))
source(paste0(scripts_dir, "optimizeGauss.R"))
source(paste0(scripts_dir, "ident.hires.noise.HPC.R"))
source(paste0(scripts_dir, "elementInfo.R"))
source(paste0(scripts_dir, "globalAssignments.HPC.R"))

# determine scan mode
if (grepl("_pos", peakgrouplist_file)) {
  scanmode <- "positive"
} else if (grepl("_neg", peakgrouplist_file)) {
  scanmode <- "negative"
}

# get replication pattern for sample names
pattern_file <- paste0(scanmode, "_repl_pattern.RData")
repl_pattern <- get(load(pattern_file))

# load peak group list and determine output file name
outpgrlist_identified <- get(load(peakgrouplist_file))

outputfile_name <- gsub(".RData", "_filled.RData", peakgrouplist_file)

# replace missing values (zeros) with random noise
peakgrouplist_filled <- replaceZeros(outpgrlist_identified, repl_pattern, scanmode, resol, outdir, thresh, ppm)

# save output
save(peakgrouplist_filled, file = outputfile_name)
