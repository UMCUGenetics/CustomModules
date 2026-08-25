# load required packages
library("ggplot2")
library("reshape2")
library("openxlsx")
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

project <- cmd_args[1]
hmdb_rlvc_file <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])
export_scripts_dir <- cmd_args[4]
path_metabolite_groups <- cmd_args[5]

# load in function scripts
source(paste0(export_scripts_dir, "generate_excel_functions.R"))

# set the number of digits for floats
options(digits = 16)

# Initialise
plot <- TRUE
export <- TRUE
control_label <- "C"
case_label <- "P"
# percentage of outliers to remove from calculation of robust scaler
perc <- 5
# Z-score for removing outliers with grubbs test
outlier_threshold <- 2

# load HMDB rlvnc table
load(hmdb_rlvc_file)

# load adduct sums object (outlist)
load("AdductSums_combined.RData")

## DrugDB output
# Filter out drugs and other metabolites with CHEMBL annotation
outlist_drugdb <- outlist[grep("CHEMBL", outlist$HMDB_ID_all), ]
if (nrow(outlist_drugdb) > 0) {
  # Add HMDB_code column with all the CHEMBL ID and an empty description column
  outlist_drugdb <- cbind(outlist_drugdb, HMDB_code = rownames(outlist_drugdb), descr = NA)
  # sort on CHEMBL ID
  outlist_drugdb <- outlist_drugdb[order(outlist_drugdb[, "HMDB_code"]), ]
  
  if (z_score == 1) {
    # calculate Z-scores with outliers removed
    outlist_drugdb_zscores <- calculate_zscores(
      outlist_drugdb, "_OutlierRemovedZscore", outlier_threshold,
      control_label, case_label
    )
    # get indices for intensity columns based on control_label and case_label
    intensity_col_ids <- get_intensity_col_index(outlist_drugdb_zscores, control_label, case_label)
    # Create Excel for drugs
    create_excel_output(outlist_drugdb_zscores, intensity_col_ids, "Drugs_", z_score, project)
  } else if (z_score == 0) {
    # get indices for intensity columns; all columns except those containing HMDB in the name
    intensity_col_ids <- 1:ncol(outlist_drugdb)
    intensity_col_ids <- intensity_col_ids[-grep("HMDB", colnames(outlist_drugdb))]
    # Create Excel for drugs
    create_excel_output(outlist_drugdb, intensity_col_ids, "Drugs_", z_score, project)
  }
}

## Filtered metabolite output
# Filter for biological relevance
peaks_in_list <- which(rownames(outlist) %in% rlvnc$HMDB_key)
outlist_subset <- outlist[peaks_in_list, ]
outlist_subset$HMDB_code <- rownames(outlist_subset)
outlist_subset$HMDB_key <- rownames(outlist_subset)
outlist_filtered <- outlist_subset %>%
  left_join(rlvnc %>% rename(sec_HMDB_ID_rlvnc = sec_HMDB_ID), by = "HMDB_key")
rownames(outlist_filtered) <- outlist_filtered$HMDB_key
# filter out all irrelevant HMDBs
outlist_filtered <- outlist_filtered %>%
  tibble::rownames_to_column("rowname") %>%
  filter(grepl("relevant|Onbekend|Internal", relevance)) %>%
  tibble::column_to_rownames("rowname")
# sort on HMDB_key
outlist_filtered <- outlist_filtered[order(outlist_filtered[, "HMDB_key"]), ]

if (z_score == 1) {
  # calculate Z-scores with outliers removed
  outlist_filtered_zscores <- calculate_zscores(
    outlist_filtered, "_OutlierRemovedZscore", outlier_threshold,
    control_label, case_label
  )
  # get indices for intensity columns
  intensity_col_ids <- get_intensity_col_index(outlist_filtered_zscores, control_label, case_label)
  control_col_idx <- get_intensity_col_index(outlist_filtered_zscores, control_label, "none")
  # Create Excel for biologically relevant metabolites
  create_excel_output(outlist_filtered_zscores, intensity_col_ids, "", z_score, project)
  # save outlist for GenerateQC step
  save(outlist_filtered_zscores, file = "outlist.RData")
  # output filtered metabolites after removal of outliers
  save_to_rdata_and_txt(outlist_filtered_zscores, "AdductSums_filtered_outliersremovedZ")
  # calculate robust Z-scores
  outlist_robust_zscore <- calculate_zscores(outlist_filtered, "_RobustZscore", perc, control_label, case_label)
  # output filtered metabolites with robust scaled Zscores
  save_to_rdata_and_txt(outlist_robust_zscore, "AdductSums_filtered_robustZ")
  # calculate Z-scores without outlier removal
  outlist <- calculate_zscores(outlist_filtered, "_Zscore", NULL, control_label, case_label)
  # output metabolites filtered on relevance
  save_to_rdata_and_txt(outlist, "AdductSums_filtered_Zscores")
} else if (z_score == 0) {
  create_excel_output(outlist_filtered, intensity_col_ids, "", z_score, project)
}

## Helix output
# get Helix IDs for extra Excel file
metabolite_files <- list.files(
  path = paste(path_metabolite_groups, "Diagnostics", sep = "/"),
  pattern = "*.txt", full.names = FALSE, recursive = FALSE
)
metab_df_helix <- NULL
for (file_index in seq_along(metabolite_files)) {
  infile <- metabolite_files[file_index]
  metab_list <- read.table(paste(path_metabolite_groups, "Diagnostics", infile, sep = "/"),
                           sep = "\t", header = TRUE, quote = ""
  )
  metab_df_helix <- rbind(metab_df_helix, metab_list)
}
# get Helix metabolites and unique HMDB IDs and remove ratio HMDBs containing A or L
metab_df_helix <- metab_df_helix %>%
  filter(Helix == "ja") %>%
  select(c(HMDB_code, HMDB_name)) %>%
  rename(H_Name = HMDB_name)
metab_list_helix <- unique(metab_df_helix$HMDB_code)
metab_list_helix <- grep("[AL]", metab_list_helix, value = TRUE, invert = TRUE)

if (z_score == 1) {
  # get intensities for Helix metabolites from dataset
  outlist_helix <- outlist_filtered_zscores %>%
    filter(HMDB_key %in% metab_list_helix) %>%
    left_join(., metab_df_helix, by = join_by(HMDB_code == HMDB_code)) %>%
    select(
      -c(HMDB_key, sec_HMDB_ID_rlvnc, name, relevance, descr, origin, fluids, tissue, disease, pathway, monositopic_mass, molecular_formula) #,
    ) 
} else if (z_score == 0) {
  # get intensities for Helix metabolites from dataset
  outlist_helix <- outlist_filtered %>%
    filter(HMDB_key %in% metab_list_helix) %>%
    left_join(., metab_df_helix, by = join_by(HMDB_code == HMDB_code)) %>%
    select(
      -c(HMDB_key, sec_HMDB_ID_rlvnc, name, relevance, descr, origin, fluids, tissue, disease, pathway, monositopic_mass, molecular_formula) #,
    ) 
}
# Create Excel for Helix
create_excel_output(outlist_helix, intensity_col_ids, "Helix_", z_score, project)

