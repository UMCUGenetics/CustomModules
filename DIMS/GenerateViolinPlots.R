# load packages
suppressPackageStartupMessages(library("dplyr"))
library(reshape2)
library(openxlsx)
library(ggplot2)
suppressPackageStartupMessages(library("gridExtra"))
library(stringr)

options(digits = 16)

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

run_name <- cmd_args[1]
export_scripts_dir <- cmd_args[2]
path_metabolite_groups <- cmd_args[3]
file_ratios_metabolites <- cmd_args[4]
file_expected_biomarkers_iem <- cmd_args[5]
file_explanation <- cmd_args[6]

# load functions
source(paste0(export_scripts_dir, "generate_violin_plots_functions.R"))
# load dataframe with intensities and Z-scores for all samples
intensities_zscore_df <- get(load("outlist.RData"))
rm(outlist)
# read input files
metabolites_ratios_df <- read.csv(file_ratios_metabolites, sep = ";", stringsAsFactors = FALSE)
expected_biomarkers_df <- read.csv(file_expected_biomarkers_iem, sep = ";", stringsAsFactors = FALSE)
expected_biomarkers_df <- expected_biomarkers_df %>%
  rename(HMDB_code = HMDB.code,
         HMDB_name = Metabolite)
explanation_violin_plot <- readLines(file_explanation)

# Set global variables
iem_variables <- list(
  top_number_iem_diseases = 5,
  threshold_iem = 5
)
top_number_iem_diseases <- 5 # number of diseases that score highest in algorithm to plot
threshold_iem <- 5 # probability score cut-off for plotting the top diseases
nr_plots_perpage <- 20 # number of violin plots per page in PDF
zscore_cutoff <- 5
protocol_name <- "DIMS_PL_DIAG"
number_of_metabolites <- list(
  highest = 20,
  lowest = 10
)

control_ids <- get_colnames_samples(intensities_zscore_df, "C")
patient_ids <- get_colnames_samples(intensities_zscore_df, "P")
all_sample_ids <- c(control_ids, patient_ids)
number_of_samples <- list(
  controls = length(control_ids),
  patients = length(patient_ids)
)

# Add Z-scores for ratios to intensities_zscore_df dataframe
intensities_zscore_ratios_df <- add_zscores_ratios_to_df(intensities_zscore_df, metabolites_ratios_df, all_sample_ids)
# for debugging:
save(intensities_zscore_ratios_df, file = "./outlist_with_ratios.RData")

# Select only the cols with zscores of the patients
zscore_patients_df <- intensities_zscore_ratios_df %>%
  select(HMDB_code, HMDB_name, any_of(paste0(patient_ids, "_Zscore"))) %>%
  rename_with(~ str_remove(.x, "_Zscore"), .cols = contains("_Zscore"))
zscore_controls_df <- intensities_zscore_ratios_df %>%
  select(HMDB_code, HMDB_name, any_of(paste0(control_ids, "_Zscore"))) %>%
  rename_with(~ str_remove(.x, "_Zscore"), .cols = contains("_Zscore"))

#### Make violin plots #####
make_and_save_violin_plot_pdfs(
  zscore_patients_df,
  zscore_controls_df,
  path_metabolite_groups,
  nr_plots_perpage,
  number_of_samples,
  run_name,
  protocol_name,
  explanation_violin_plot,
  number_of_metabolites
)

#### Run the IEM algorithm #########
diem_probability_score <- run_diem_algorithm(expected_biomarkers_df, zscore_patients_df, patient_ids)

save_prob_scores_to_excel(diem_probability_score, run_name)

#### Generate dIEM plots #########
patient_no_iem <- make_and_save_diem_plots(
  diem_probability_score,
  patient_ids,
  expected_biomarkers_df,
  zscore_patients_df,
  zscore_controls_df,
  nr_plots_perpage,
  number_of_samples,
  number_of_metabolites,
  iem_variables
)

if (length(patient_no_iem) > 0) {
  save_patient_no_iem(iem_variables$threshold_iem, patient_no_iem)
}
