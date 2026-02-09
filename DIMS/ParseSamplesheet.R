# define parameters
args <- commandArgs(trailingOnly = TRUE)

sample_sheet <- as.data.frame(read.csv(args[1], sep = "\t"))
preprocessing_scripts_dir <- args[2]

# load in function script
source(paste0(preprocessing_scripts_dir, "parse_samplesheet_functions.R"))

# generate the replication pattern
repl_pattern <- generate_repl_pattern(sample_sheet)

# write the replication pattern to text file for troubleshooting purposes
sink("replication_pattern.txt")
print(repl_pattern)
sink()

# save replication pattern to file
save(repl_pattern, file = "init.RData")
