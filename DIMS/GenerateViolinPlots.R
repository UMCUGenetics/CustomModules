# load packages
suppressPackageStartupMessages(library("dplyr"))
library(reshape2)
library(openxlsx)
library(ggplot2)
suppressPackageStartupMessages(library("gridExtra"))
library(stringr)

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
# read input files
ratios_metabs_df <- read.csv(file_ratios_metabolites, sep = ";", stringsAsFactors = FALSE)
expected_biomarkers_df <- read.csv(file_expected_biomarkers_iem, sep = ";", stringsAsFactors = FALSE)
explanation_violin_plot <- readLines(file_explanation)


## Set global variables
output_dir <- "./"      # path: output folder for dIEM and violin plots
top_number_iem_diseases <- 5         # number of diseases that score highest in algorithm to plot
threshold_iem <- 5      # probability score cut-off for plotting the top diseases
ratios_cutoff <- -5     # z-score cutoff of axis on the left for top diseases
nr_plots_perpage <- 20  # number of violin plots per page in PDF
zscore_cutoff <- 5
xaxis_cutoff <- 20
protocol_name <- "DIMS_PL_DIAG"

# Remove columns, move HMDB_code & HMDB_name column to the front, change intensity columns to numeric
intensities_zscore_df <- intensities_zscore_df %>%
  select(-c(plots, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID, HMDB_key, sec_HMBD_ID_rlvnc, name, 
            relevance, descr, origin, fluids, tissue, disease, pathway, nr_ctrls)) %>%
  relocate(c(HMDB_code, HMDB_name)) %>%
  rename(mean_controls = avg_ctrls, sd_controls = sd_ctrls) %>%
  mutate(across(!c(HMDB_name, HMDB_code), as.numeric))

# Get the controls and patient IDs, select the intensity columns
controls <- colnames(intensities_zscore_df)[grepl("^C", colnames(intensities_zscore_df)) &
                                              !grepl("_Zscore$", colnames(intensities_zscore_df))]
control_intensities_cols_index <- which(colnames(intensities_zscore_df) %in% controls)
nr_of_controls <- length(controls)

patients <- colnames(intensities_zscore_df)[grepl("^P", colnames(intensities_zscore_df)) &
                                              !grepl("_Zscore$", colnames(intensities_zscore_df))]
patient_intensities_cols_index <- which(colnames(intensities_zscore_df) %in% patients)
nr_of_patients <- length(patients)

intensity_cols_index <- c(control_intensities_cols_index, patient_intensities_cols_index)
intensity_cols <- colnames(intensities_zscore_df)[intensity_cols_index]

#### Calculate ratios of intensities for metabolites ####
# Prepare empty data frame to fill with ratios
ratio_zscore_df <- data.frame(matrix(
  ncol = ncol(intensities_zscore_df),
  nrow = nrow(ratios_metabs_df)
))
colnames(ratio_zscore_df) <- colnames(intensities_zscore_df)

# put HMDB info into first two columns of ratio_zscore_df
ratio_zscore_df$HMDB_code <- ratios_metabs_df$HMDB.code
ratio_zscore_df$HMDB_name <- ratios_metabs_df$Ratio_name

for (row_index in 1:nrow(ratios_metabs_df)) {
  numerator_intensities <- get_intentities_for_ratios(ratios_metabs_df, row_index, 
                                                  intensities_zscore_df, "HMDB_numerator", intensity_cols)
  denominator_intensities <- get_intentities_for_ratios(ratios_metabs_df, row_index, 
                                                    intensities_zscore_df, "HMDB_denominator", intensity_cols)
  # calculate intensity ratios
  ratio_zscore_df[row_index, intensity_cols_index] <- log2(numerator_intensities / denominator_intensities)
}
# Calculate means and SD's of the calculated ratios for Controls
ratio_zscore_df[, "mean_controls"] <- apply(ratio_zscore_df[, control_intensities_cols_index], 1, mean)
ratio_zscore_df[, "sd_controls"]   <- apply(ratio_zscore_df[, control_intensities_cols_index], 1, sd)

# Calculate Zscores for the ratios
samples_zscore_columns <- get_zscore_columns(colnames(intensities_zscore_df), intensity_cols)
ratio_zscore_df[, samples_zscore_columns] <- (ratio_zscore_df[, intensity_cols] - ratio_zscore_df[, "mean_controls"]) /
  ratio_zscore_df[, "sd_controls"]

intensities_zscore_ratios_df <- rbind(intensities_zscore_df, ratio_zscore_df)

# for debugging:
save(intensities_zscore_ratios_df, file = paste0(output_dir, "/outlist_with_ratios.RData"))

# Select only the cols with zscores of the patients
zscore_patients_df <- intensities_zscore_ratios_df %>% select(HMDB_code, HMDB_name, any_of(paste0(patients, "_Zscore")))
zscore_controls_df <- intensities_zscore_ratios_df %>% select(HMDB_code, HMDB_name, any_of(paste0(controls, "_Zscore")))

#### Make violin plots #####
# preparation
colnames(zscore_patients_df) <- gsub("_Zscore", "", colnames(zscore_patients_df))
colnames(zscore_controls_df) <- gsub("_Zscore", "", colnames(zscore_controls_df))

expected_biomarkers_df <- expected_biomarkers_df %>% rename(HMDB_code = HMDB.code, HMDB_name = Metabolite)

expected_biomarkers_info <- expected_biomarkers_df %>% 
  select(c(Disease, HMDB_code, HMDB_name)) %>% 
  distinct(Disease, HMDB_code, .keep_all = TRUE)

metabolite_dirs <- list.files(path = path_metabolite_groups, full.names = FALSE, recursive = FALSE)
for (metabolite_dir in metabolite_dirs) {
  # create a directory for the output PDFs
  pdf_dir <- paste(output_dir, metabolite_dir, sep = "/")
  dir.create(pdf_dir, showWarnings = FALSE)

  metab_list_all <- get_list_metabolites(paste(path_metabolite_groups, metabolite_dir, sep = "/"))

  # prepare list of metabolites; max nr_plots_perpage on one page
  metab_interest_sorted <- combine_metab_info_zscores(metab_list_all, zscore_patients_df)
  metab_interest_controls <- combine_metab_info_zscores(metab_list_all, zscore_controls_df)
  metab_perpage <- prepare_data_perpage(metab_interest_sorted, metab_interest_controls,
                                        nr_plots_perpage, nr_of_patients, nr_of_controls)

  # for Diagnostics metabolites to be saved in Helix
  if(grepl("Diagnost", pdf_dir)) {
    # get table that combines DIMS results with stofgroepen/Helix table
    dims_helix_table <- get_patient_data_to_helix(metab_interest_sorted, metab_list_all)
    
    # check if run contains Diagnostics patients (e.g. "P2024M"), not for research runs
    if(any(is_diagnostic_patient(dims_helix_table$Sample))){
      # get output file for Helix
      output_helix <- output_for_helix(protocol_name, dims_helix_table)
      # write output to file
      path_helixfile <- paste0(output_dir, "output_Helix_", run_name,".csv")
      write.csv(output_helix, path_helixfile, quote = F, row.names = F)
    }
  }
  
  # make violin plots per patient
  for (patient_id in patients) {
    # for category Diagnostics, make list of metabolites that exceed alarm values for this patient
    # for category Other, make list of top highest and lowest Z-scores for this patient
    if (grepl("Diagnost", pdf_dir)) {
      top_metabs_patient <- prepare_alarmvalues(patient_id, dims_helix_table)
    } else {
      top_metabs_patient <- prepare_toplist(patient_id, zscore_patients_df)
    }

    # generate normal violin plots
    create_pdf_violin_plots(pdf_dir, patient_id, metab_perpage, top_metabs_patient, explanation_violin_plot)
  }

}

#### Run the IEM algorithm #########
expected_biomarkers_df <- expected_biomarkers_df %>% rename(HMDB_code = HMDB.code, HMDB_name = Metabolite)

diem_probability_score <- run_diem_algorithm(expected_biomarkers_df, zscore_patients_df, patients)

save_prob_scores_to_Excel(diem_probability_score, output_dir, run_name)


#### Generate dIEM plots #########
diem_plot_dir <- paste(output_dir, "dIEM_plots", sep = "/")
dir.create(diem_plot_dir)

colnames(diem_probability_score) <- gsub("_Zscore", "", colnames(diem_probability_score))
patient_no_iem <- c()

for (patient_id in patients) {
  # Select the top IEMs and filter on the IEM threshold
  patient_top_iems_probs <- diem_probability_score %>%
    select(c(Disease, !!sym(patient_id))) %>%
    arrange(desc(!!sym(patient_id))) %>%
    slice(1:top_number_iem_diseases) %>%
    filter(!!sym(patient_id) >= threshold_iem)
  
  if (nrow(patient_top_iems_probs) > 0) {
    top_iems <- patient_top_iems_probs %>% pull(Disease)
    # Get the metabolites for each IEM and their probability
    metabs_iems_names <- c()
    metabs_iems <- lapply(top_iems, function(iem) {
      iem_probablity <- patient_top_iems_probs %>% filter(Disease == iem) %>% pull(!!sym(patient_id))
      metabs_iems_names <- c(metabs_iems_names, paste0(iem, ", probability score ", iem_probablity))
      metab_iem <- expected_biomarkers_df %>% filter(Disease == iem) %>% select(HMDB_code, HMDB_name)
      return(metab_iem)
    })
    names(metabs_iems) <- metabs_iems_names
    
    # Get the Z-scores with metabolite information
    metab_iem_sorted <- combine_metab_info_zscores(metabs_iems, zscore_patients_df)
    metab_iem_controls <- combine_metab_info_zscores(metabs_iems, zscore_controls_df)
    # Get a list of dataframes for each IEM
    diem_metab_perpage <- prepare_data_perpage(metab_iem_sorted, metab_iem_controls,
                                               nr_plots_perpage, nr_of_patients, nr_of_controls)
    # Get a dataframe of the top metabolites 
    top_metabs_patient <- prepare_toplist(patient_id, zscore_patients_df)
    
    # Generate and save dIEM violin plots
    create_pdf_violin_plots(diem_plot_dir, patient_id, diem_metab_perpage, top_metabs_patient, explanation_violin_plot)
    
  } else {
    patient_no_iem <- c(patient_no_iem, patient_id)
  }
}

if (length(patient_no_iem) > 0) {
  patient_no_iem <- c(paste0("The following patient(s) did not have dIEM probability scores higher than ", threshold_iem, " :"),
                      patient_no_iem)
  write(file = paste0(output_dir, "missing_probability_scores.txt"), patient_no_iem)
}
