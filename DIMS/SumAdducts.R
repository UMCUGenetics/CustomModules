# load libraries
library("dplyr")
library("argparse")

parser <- ArgumentParser(description = "SumAdducts")

parser$add_argument("--hmdbpart_main_file", dest = "hmdbpart_main_file",
                    help = "RData file containing part of HMDB database (only main metabolites)", required = TRUE)
parser$add_argument("--preprocessing_scripts_dir", dest = "preprocessing_scripts_dir",
                    help = "File path to the directory containing functions used", required = TRUE)
parser$add_argument("--z_score", dest = "z_score",
                    help = "Boolean indicating whether Z-scores should be calculated (1) or not (0)", required = TRUE)

args <- parser$parse_args()

# define parameters
hmdbpart_main_file <- args$hmdbpart_main_file
z_score <- as.numeric(args$z_score)

# load in function scripts
source(paste0(args$preprocessing_scripts_dir, "sum_intensities_adducts.R"))

if (grepl("positive_hmdb", hmdbpart_main_file)) {
  scanmode <- "positive"
  # for the adduct sum: include adducts M+Na (1) and M+K (2)
  adducts <- c(1, 2)
} else if (grepl("negative_hmdb", hmdbpart_main_file)) {
  scanmode <- "negative"
  # for the adduct sum: include adduct M+Cl (1)
  adducts <- c(1)
}

# load input files
collect_file <- paste0("outlist_identified_", scanmode, ".RData")
peakgroup_list <- get(load(collect_file))
hmdb_main_part <- get(load(hmdbpart_main_file))

# get the number from the file name
batch_number <- strsplit(basename(hmdbpart_main_file), ".", fixed = TRUE)[[1]][2]

# sum adducts and save output
adductsum <- sum_intensities_adducts(peakgroup_list, hmdb_main_part, adducts, z_score)
save(adductsum, file = paste(scanmode, "_", batch_number, "_SummedAdducts.RData", sep = ""))
