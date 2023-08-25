#!/usr/bin/Rscript
# adapted from 9-runFillMissing.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

# define parameters 
peakgrouplist_file <- cmd_args[1]
#pattern_file       <- cmd_args[2]
scripts_dir        <- cmd_args[2]
#thresh  <- as.numeric(cmd_args[4])
#resol   <- as.numeric(cmd_args[5])
#ppm     <- as.numeric(cmd_args[6])
thresh  <- cmd_args[3]
resol   <- cmd_args[4]
ppm     <- cmd_args[5]
# outdir <- cmd_args[2]
# scanmode <- cmd_args[3]
outdir <- "./"
if (grepl("_pos", peakgrouplist_file)) { scanmode = "positive" } else
    if (grepl("_neg", peakgrouplist_file)) { scanmode = "negative" }

# load in function scripts
source(paste0(scripts_dir, "AddOnFunctions/replaceZeros.R"))
source(paste0(scripts_dir, "AddOnFunctions/generateGaussian.R"))
source(paste0(scripts_dir, "AddOnFunctions/ident.hires.noise.HPC.R"))

# for debugging:
print(peakgrouplist_file)
#print(pattern_file)
print(scripts_dir)
print(thresh)
print(resol)
print(ppm)

outputfile_name <- gsub(".RData", "_filled.RData", peakgrouplist_file)
print(outputfile_name)

# get replication pattern for sample names
pattern_file <- paste0(scanmode, "_repl_pattern.RData")
repl_pattern <- get(load(pattern_file))
print(repl_pattern[[1]])

# peakgrouplist <- get(load(peakgrouplist_file))
load(peakgrouplist_file)
batch_number <- strsplit(basename(peakgrouplist_file), ".",fixed = TRUE)[[1]][1]

print(dim(outpgrlist.identified))
print(colnames(outpgrlist.identified))
print(batch_number)

# replace missing values (zeros) with random noise
peakgrouplist_filled <- replaceZeros(outpgrlist.identified, repl_pattern, scanmode, resol, outdir, thresh, ppm)

# save output
save(peakgrouplist_filled, file=paste0("./", outputfile_name)
save(peakgrouplist_filled, file=paste0("./", batch_number, scanmode, "identified_filled.RData")
