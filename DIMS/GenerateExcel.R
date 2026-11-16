# Load required packages
library("ggplot2")
library("reshape2")
library("openxlsx")
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))

# Define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

project <- cmd_args[1]
hmdb_rlvc_file <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])
export_scripts_dir <- cmd_args[4]
path_metabolite_groups <- cmd_args[5]

# Load in function scripts
source(paste0(export_scripts_dir, "/generate_excel_functions.R"))

# Initialize
options(digits = 16)
plot <- TRUE
export <- TRUE
control_label <- "C"
case_label <- "P"
# Percentage of outliers to remove from calculation of robust scaler
perc <- 5
# Z-score for removing outliers with grubbs test
outlier_threshold <- 2

# Load HMDB rlvnc table
load(hmdb_rlvc_file)

# Load peakgroup_list object
peakgroup_list <- get(load("AdductSums_combined.RData"))

## DrugDB output
# Filter out drugs and other metabolites with CHEMBL annotation
peakgroup_list_drugdb <- peakgroup_list[grep("CHEMBL", peakgroup_list$HMDB_ID_all), ]
if (nrow(peakgroup_list_drugdb) > 0) {
  # Add HMDB_code column with all the CHEMBL ID and an empty description column
  peakgroup_list_drugdb <- cbind(peakgroup_list_drugdb, HMDB_code = rownames(peakgroup_list_drugdb), descr = NA)
  # Sort on CHEMBL ID
  peakgroup_list_drugdb <- peakgroup_list_drugdb[order(peakgroup_list_drugdb[, "HMDB_code"]), ]
  
  if (z_score == 1) {
    # Calculate Z-scores with outliers removed
    peakgroup_list_drugdb_zscores <- calculate_zscores(
      peakgroup_list_drugdb, "_OutlierRemovedZscore", outlier_threshold,
      control_label, case_label
    )
    # Get indices for intensity columns based on control_label and case_label
    intensity_col_ids <- get_intensity_col_index(peakgroup_list_drugdb_zscores, control_label, case_label)
    # Create Excel for drugs
    create_excel_output(peakgroup_list_drugdb_zscores, intensity_col_ids, "Drugs_", z_score, project)
    # save to RData file
    save(peakgroup_list_drugdb_zscores, file = "Drugs_in_dataset.RData")
  } else if (z_score == 0) {
    # Get indices for intensity columns; all columns except those containing HMDB in the name
    intensity_col_ids <- 1:ncol(peakgroup_list_drugdb)
    intensity_col_ids <- intensity_col_ids[-grep("HMDB", colnames(peakgroup_list_drugdb))]
    # Create Excel for drugs
    create_excel_output(peakgroup_list_drugdb, intensity_col_ids, "Drugs_", z_score, project)
  }
}

## Filtered metabolite output
# Filter for biological relevance
peaks_in_list <- which(rownames(peakgroup_list) %in% rlvnc$HMDB_key)
peakgroup_list_subset <- peakgroup_list[peaks_in_list, ]
peakgroup_list_subset$HMDB_code <- rownames(peakgroup_list_subset)
peakgroup_list_subset$HMDB_key <- rownames(peakgroup_list_subset)
peakgroup_list_filtered <- peakgroup_list_subset %>%
  left_join(rlvnc %>% rename(sec_HMDB_ID_rlvnc = sec_HMDB_ID), by = "HMDB_key")
rownames(peakgroup_list_filtered) <- peakgroup_list_filtered$HMDB_key
# Filter out all irrelevant HMDBs
peakgroup_list_filtered <- peakgroup_list_filtered %>%
  tibble::rownames_to_column("rowname") %>%
  filter(grepl("relevant|Onbekend|Internal", relevance)) %>%
  tibble::column_to_rownames("rowname")
# Sort on HMDB_key
peakgroup_list_filtered <- peakgroup_list_filtered[order(peakgroup_list_filtered[, "HMDB_key"]), ]

if (z_score == 1) {
  # Calculate Z-scores with outliers removed
  peakgroup_list_filtered_zscores <- calculate_zscores(
    peakgroup_list_filtered, "_OutlierRemovedZscore", outlier_threshold,
    control_label, case_label
  )
  # Get indices for intensity columns
  intensity_col_ids <- get_intensity_col_index(peakgroup_list_filtered_zscores, control_label, case_label)
  control_col_idx <- get_intensity_col_index(peakgroup_list_filtered_zscores, control_label, "none")
  # Create Excel for biologically relevant metabolites
  create_excel_output(peakgroup_list_filtered_zscores, intensity_col_ids, "", z_score, project)
  # Save peakgroup_list for GenerateQC step
  save(peakgroup_list_filtered_zscores, file = "peakgrouplist.RData")
  # Output filtered metabolites after removal of outliers
  save_to_rdata_and_txt(peakgroup_list_filtered_zscores, "AdductSums_filtered_outliersremovedZ")
  # Calculate robust Z-scores
  peakgroup_list_robust_zscores <- calculate_zscores(peakgroup_list_filtered, "_RobustZscore", perc, control_label, case_label)
  # Output filtered metabolites with robust scaled Zscores
  save_to_rdata_and_txt(peakgroup_list_robust_zscores, "AdductSums_filtered_robustZ")
  # Calculate Z-scores without outlier removal
  peakgroup_list_zscores <- calculate_zscores(peakgroup_list_filtered, "_Zscore", NULL, control_label, case_label)
  # Output metabolites filtered on relevance
  save_to_rdata_and_txt(peakgroup_list_zscores, "AdductSums_filtered_Zscores")
} else if (z_score == 0) {
  create_excel_output(peakgroup_list_filtered, intensity_col_ids, "", z_score, project)
}

## Helix output
# Get Helix IDs for extra Excel file
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
# Get Helix metabolites and unique HMDB IDs and remove ratio HMDBs containing A or L
metab_df_helix <- metab_df_helix %>%
  filter(Helix == "ja") %>%
  select(c(HMDB_code, HMDB_name)) %>%
  rename(H_Name = HMDB_name)
metab_list_helix <- unique(metab_df_helix$HMDB_code)
metab_list_helix <- grep("[AL]", metab_list_helix, value = TRUE, invert = TRUE)

if (z_score == 1) {
  # Get intensities for Helix metabolites from dataset
  peakgroup_list_helix <- peakgroup_list_filtered_zscores %>%
    filter(HMDB_key %in% metab_list_helix) %>%
    left_join(., metab_df_helix, by = join_by(HMDB_code == HMDB_code)) %>%
    select(
      -c(HMDB_key, sec_HMDB_ID_rlvnc, name, relevance, descr, origin, fluids, tissue, disease, pathway, monositopic_mass, molecular_formula) #,
    ) 
} else if (z_score == 0) {
  # Get intensities for Helix metabolites from dataset
  peakgroup_list_helix <- peakgroup_list_filtered %>%
    filter(HMDB_key %in% metab_list_helix) %>%
    left_join(., metab_df_helix, by = join_by(HMDB_code == HMDB_code)) %>%
    select(
      -c(HMDB_key, sec_HMDB_ID_rlvnc, name, relevance, descr, origin, fluids, tissue, disease, pathway, monositopic_mass, molecular_formula) #,
    ) 
}

# Create Excel for Helix
create_excel_output(peakgroup_list_helix, intensity_col_ids, "Helix_", z_score, project)
