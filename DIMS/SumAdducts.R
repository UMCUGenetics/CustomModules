library("dplyr")
# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

hmdbpart_main_file <- cmd_args[1]
scripts_dir <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])

# load in function scripts
source(paste0(scripts_dir, "sum_intensities_adducts.R"))

if (grepl("positive_hmdb", hmdbpart_main_file)) {
  scanmode <- "positive"
} else if (grepl("negative_hmdb", hmdbpart_main_file)) {
  scanmode <- "negative"
}

# load input files
collect_file <- paste0("outlist_identified_", scanmode, ".RData")
peakgroup_list <- get(load(collect_file))
hmdb_main_part <- get(load(hmdbpart_main_file))

# get the number from the file name
batch_number <- strsplit(basename(hmdbpart_main_file), ".", fixed = TRUE)[[1]][2]

# sum adducts and save output
adductsum <- sum_intensities_adducts_perclass(peakgroup_list, hmdb_main_part, z_score)
save(adductsum, file = paste(scanmode, "_", batch_number, "_SummedAdducts.RData", sep = ""))
