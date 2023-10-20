#!/usr/bin/Rscript
# adapted from 9-runFillMissing.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

# define parameters 
peakgrouplist_file <- cmd_args[1]
scripts_dir        <- cmd_args[2]
thresh  <- as.numeric(cmd_args[3])
resol   <- as.numeric(cmd_args[4])
ppm     <- as.numeric(cmd_args[5])
outdir <- "./"

print(peakgrouplist_file)

if (grepl("_pos", peakgrouplist_file)) { scanmode = "positive" } else
    if (grepl("_neg", peakgrouplist_file)) { scanmode = "negative" }

print(scanmode)
print(scripts_dir)

# load in function scripts
source(paste0(scripts_dir, "AddOnFunctions/replaceZeros.R"))
source(paste0(scripts_dir, "AddOnFunctions/generateGaussian.R"))
source(paste0(scripts_dir, "AddOnFunctions/getFwhm.R"))
source(paste0(scripts_dir, "AddOnFunctions/getSD.R"))
source(paste0(scripts_dir, "AddOnFunctions/getArea.R"))
source(paste0(scripts_dir, "AddOnFunctions/optimizeGauss.R"))
source(paste0(scripts_dir, "AddOnFunctions/ident.hires.noise.HPC.R"))
source(paste0(scripts_dir, "AddOnFunctions/elementInfo.R"))
source(paste0(scripts_dir, "AddOnFunctions/globalAssignments.HPC.R"))

# get replication pattern for sample names
pattern_file <- paste0(scanmode, "_repl_pattern.RData")
repl_pattern <- get(load(pattern_file))

print(head(repl_pattern))

# load peak group list and determine output file name
outpgrlist_identified <- get(load(peakgrouplist_file))

print(head(outpgrlist_identified))

outputfile_name <- gsub(".RData", "_filled.RData", peakgrouplist_file)

print(outputfile_name)

# replace missing values (zeros) with random noise
peakgrouplist_filled <- replaceZeros(outpgrlist_identified, repl_pattern, scanmode, resol, outdir, thresh, ppm)

# save output 
save(peakgrouplist_filled, file=paste0("./", outputfile_name))
