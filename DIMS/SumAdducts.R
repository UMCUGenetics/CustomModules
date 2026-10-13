# Load required libraries
library("dplyr")

# Define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

hmdbpart_main_file <- cmd_args[1]
scripts_dir <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])

# Load in function scripts
source(paste0(scripts_dir, "sum_intensities_adducts.R"))

if (grepl("positive_hmdb", hmdbpart_main_file)) {
  scanmode <- "positive"
  # For the adduct sum: include adducts M+Na (1) and M+K (2)
  adducts <- c(1, 2)
} else if (grepl("negative_hmdb", hmdbpart_main_file)) {
  scanmode <- "negative"
  # For the adduct sum: include adduct M+Cl (1)
  adducts <- c(1)
}

# Load input files
collect_file <- paste0("peakgroup_list_identified_", scanmode, ".RData")
peakgroup_list <- get(load(collect_file))
hmdb_main_part <- get(load(hmdbpart_main_file))

# Get the number from the file name
batch_number <- strsplit(basename(hmdbpart_main_file), ".", fixed = TRUE)[[1]][2]

# Sum adducts and save output
adductsum <- sum_intensities_adducts(peakgroup_list, hmdb_main_part, adducts, z_score)
save(adductsum, file = paste(scanmode, "_", batch_number, "_SummedAdducts.RData", sep = ""))
