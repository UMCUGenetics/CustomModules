# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

hmdbpart_main_file <- cmd_args[1]
scripts_dir <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])

# load in function scripts
source(paste0(scripts_dir, "sum_adducts.R"))

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
load(collect_file)
repl_file <- paste0(scanmode, "_repl_pattern.RData")
load(repl_file)
hmdb_main_part <- get(load(hmdbpart_main_file))

# get the number from the file name
batch_number <- strsplit(basename(hmdbpart_main_file), ".", fixed = TRUE)[[1]][2]

# sum adducts and save output
summed_list <- sum_adducts(outlist_total, hmdb_main_part, names(repl_pattern_filtered), adducts, z_score)
save(summed_list, file = paste(scanmode, "_", batch_number, "_SummedAdducts.RData", sep = ""))
