## Combining all AdductSums part files for each scanmode and
# combine intensities if present in both scanmodes
library("argparse")
suppressMessages(library("dplyr"))

parser <- ArgumentParser(description = "GenerateExcel")

parser$add_argument("--preprocessing_scripts_dir", dest = "preprocessing_scripts_dir",
                    help = "File path to the directory containing functions used in this script", required = TRUE)

args <- parser$parse_args()

source(paste0(args$preprocessing_scripts_dir, "collect_sum_adducts_functions.R"))

# collect all AdductSums part files for each scanmode and save to RData file
outlist_tot_pos <- combine_sum_adduct_parts("positive")
save(outlist_tot_pos, file = "AdductSums_positive.RData")
outlist_tot_neg <- combine_sum_adduct_parts("negative")
save(outlist_tot_neg, file = "AdductSums_negative.RData")

# combine intensities of both scanmodi
outlist <- combine_scanmodes_intensities(outlist_tot_pos, outlist_tot_neg)
save(outlist, file = "AdductSums_combined.RData")
