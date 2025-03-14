# For untargeted metabolomics, this tool calculates probability scores for
# metabolic disorders. In addition, it provides visual support with violin plots
# of the DIMS measurements for the lab specialists.
# Input needed:
# 1. Excel file in which metabolites are listed with their intensities for
#    controls (with C in samplename) and patients (with P in samplename) and their
#    corresponding Z-scores.
# 2. All files from github: https://github.com/UMCUGenetics/DIMS

## adapted from 15-dIEM_violin.R

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
scripts_dir <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])
path_metabolite_groups <- cmd_args[4]
file_ratios_metabolites <- cmd_args[5]
file_expected_biomarkers_iem <- cmd_args[6]
file_explanation <- cmd_args[7]
file_isomers <- cmd_args[8]

if (z_score == 1){
  # path: output folder for dIEM and violin plots
  output_dir <- "./"

  file.copy(file_isomers, output_dir)

  # load functions
  source(paste0(scripts_dir, "check_same_samplename.R"))
  source(paste0(scripts_dir, "prepare_data.R"))
  source(paste0(scripts_dir, "prepare_data_perpage.R"))
  source(paste0(scripts_dir, "prepare_toplist.R"))
  source(paste0(scripts_dir, "create_violin_plots.R"))
  source(paste0(scripts_dir, "prepare_alarmvalues.R"))
  source(paste0(scripts_dir, "output_helix.R"))
  source(paste0(scripts_dir, "get_patient_data_to_helix.R"))
  source(paste0(scripts_dir, "add_lab_id_and_onderzoeksnummer.R"))
  source(paste0(scripts_dir, "is_diagnostic_patient.R"))

  # number of diseases that score highest in algorithm to plot
  top_nr_iem <- 5
  # probability score cut-off for plotting the top diseases
  threshold_iem <- 5
  # z-score cutoff of axis on the left for top diseases
  ratios_cutoff <- -5
  # number of violin plots per page in PDF
  nr_plots_perpage <- 20

  # binary variable: run function, yes(1) or no(0)
  if (z_score == 1) { 
    algorithm <- ratios <- violin <- 1 
  } else {
    algorithm <- ratios <- violin <- 0 
  }	
  # are the sample names headers on row 1 or row 2 in the DIMS excel? (default 1)
  header_row <- 1
  # column name where the data starts (default B)
  col_start <- "B"
  zscore_cutoff <- 5
  xaxis_cutoff <- 20
  protocol_name <- "DIMS_PL_DIAG"

  #### STEP 1: Preparation ####
  # in: run_name, path_dims_file, header_row ||| out: output_dir, DIMS

  # load outlist instead of excel file
  load("outlist.RData")

  # save outlist as dims_xls, will be changed during refactor 
  dims_xls <- outlist
  rm(outlist)

  #### STEP 2: Edit DIMS data #####
  # in: dims_xls ||| out: Data, nr_contr, nr_pat
  # Input: the xlsx file that comes out of the pipeline with format:
  # [plots] [C] [P] [summary columns] [C_Zscore] [P_Zscore]
  # Output: "_CSV.csv" file that is suited for the algorithm in shiny.

  # Determine the number of Contols and Patients in column names:
  nr_contr <- length(grep("C", names(dims_xls))) / 2
  nr_pat   <- length(grep("P", names(dims_xls))) / 2
  # total number of samples
  nrsamples <- nr_contr + nr_pat
  # check whether the number of intensity columns equals the number of Zscore columns
  if (nr_contr + nr_pat != length(grep("_Zscore", names(dims_xls)))) {
    cat("\n**** Error: there aren't as many intensities listed as Zscores")
  }
  cat(paste0("\n\n------------\n", nr_contr, " controls \n", nr_pat, " patients\n------------\n\n"))

  # Move the columns HMDB_code and HMDB_name to the beginning.
  hmdb_info_cols <- c(which(colnames(dims_xls) == "HMDB_code"), which(colnames(dims_xls) == "HMDB_name"))
  other_cols <- seq_along(1:ncol(dims_xls))[-hmdb_info_cols]
  dims_xls_copy <- dims_xls[, c(hmdb_info_cols, other_cols)]
  # Remove the columns from 'name' to 'pathway'
  from_col <- which(colnames(dims_xls_copy) == "name")
  to_col <- which(colnames(dims_xls_copy) == "pathway")
  dims_xls_copy <- dims_xls_copy[, -c(from_col:to_col)]
  # in case the excel had an empty "plots" column, remove it
  if ("plots" %in% colnames(dims_xls_copy)) {
    dims_xls_copy <- dims_xls_copy[, -grep("plots", colnames(dims_xls_copy))]
  }
  # Rename columns
  names(dims_xls_copy) <- gsub("avg.ctrls", "Mean_controls", names(dims_xls_copy))
  names(dims_xls_copy) <- gsub("sd.ctrls",  "SD_controls", names(dims_xls_copy))
  names(dims_xls_copy) <- gsub("HMDB_code", "HMDB.code", names(dims_xls_copy))
  names(dims_xls_copy) <- gsub("HMDB_name", "HMDB.name", names(dims_xls_copy))

  # intensity columns and mean and standard deviation of controls
  numeric_cols <- c(3:ncol(dims_xls_copy))
  # make sure all values are numeric
  dims_xls_copy[, numeric_cols] <- sapply(dims_xls_copy[, numeric_cols], as.numeric)

  if (exists("dims_xls_copy") & (length(dims_xls_copy) < length(dims_xls))) {
    cat("\n### Step 2 # Edit dims data is done.\n")
  } else {
    cat("\n**** Error: Could not execute step 2 \n")
  }

  #### STEP 3: Calculate ratios of intensities for metabolites ####
  # in: ratios, file_ratios_metabolites, dims_xls_copy, nr_contr, nr_pat ||| out: Zscore (+file)
  # This script loads the file with Ratios (file_ratios_metabolites) and calculates
  # the ratios of the intensities of the given metabolites. It also calculates
  # Zs-cores based on the avg and sd of the ratios of the controls.

  # Input: dataframe with intenstities and Zscores of controls and patients:
  # [HMDB.code] [HMDB.name] [C] [P] [Mean_controls] [SD_controls] [C_Zscore] [P_Zscore]

  # Output: "_CSV.csv" file that is suited for the algorithm, with format:
  # "_Ratios_CSV.csv" file, same file as above, but with ratio rows added.

  if (ratios == 1) {
    cat(paste0("\nloading ratios file:\n ->  ", file_ratios_metabolites, "\n"))
    ratio_input <- read.csv(file_ratios_metabolites, sep = ";", stringsAsFactors = FALSE)

    # Prepare empty data frame to fill with ratios
    ratio_list <- setNames(data.frame(matrix(
      ncol = ncol(dims_xls_copy),
      nrow = nrow(ratio_input)
    )), colnames(dims_xls_copy))
    ratio_list <- as.data.frame(ratio_list)

    # put HMDB info into first two columns of ratio_list
    ratio_list[, 1:2] <- ratio_input[, 1:2]

    # look for intensity columns (exclude Zscore columns)
    control_cols <- grep("C", colnames(ratio_list)[1:which(colnames(ratio_list) == "Mean_controls")])
    patient_cols <- grep("P", colnames(ratio_list)[1:which(colnames(ratio_list) == "Mean_controls")])
    intensity_cols <- c(control_cols, patient_cols)
    # calculate each of the ratios of intensities
    for (ratio_index in 1:nrow(ratio_input)) {
      ratio_numerator <- ratio_input[ratio_index, "HMDB_numerator"]
      ratio_numerator <- strsplit(ratio_numerator, "plus")[[1]]
      ratio_denominator <- ratio_input[ratio_index, "HMDB_denominator"]
      ratio_denominator <- strsplit(ratio_denominator, "plus")[[1]]
      # find these HMDB IDs in dataset. Could be a sum of multiple metabolites
      sel_denominator <- sel_numerator <- c()
      for (numerator_index in 1:length(ratio_numerator)) {
        sel_numerator <- c(sel_numerator, which(dims_xls_copy[, "HMDB.code"] == ratio_numerator[numerator_index]))
      }
      for (denominator_index in 1:length(ratio_denominator)) {
        # special case for sum of metabolites (dividing by one)
        if (ratio_denominator[denominator_index] != "one") {
          sel_denominator <- c(sel_denominator, which(dims_xls_copy[, "HMDB.code"] == ratio_denominator[denominator_index]))
        }
      }
      # calculate ratio
      if (ratio_denominator[denominator_index] != "one") {
        ratio_list[ratio_index, intensity_cols] <- apply(dims_xls_copy[sel_numerator, intensity_cols], 2, sum) /
          apply(dims_xls_copy[sel_denominator, intensity_cols], 2, sum)
      } else {
        # special case for sum of metabolites (dividing by one)
        ratio_list[ratio_index, intensity_cols] <- apply(dims_xls_copy[sel_numerator, intensity_cols], 2, sum)
      }
      # calculate log of ratio
      ratio_list[ratio_index, intensity_cols] <- log2(ratio_list[ratio_index, intensity_cols])
    }

    # Calculate means and SD's of the calculated ratios for Controls
    ratio_list[, "Mean_controls"] <- apply(ratio_list[, control_cols], 1, mean)
    ratio_list[, "SD_controls"]   <- apply(ratio_list[, control_cols], 1, sd)

    # Calc z-scores with the means and SD's of Controls
    zscore_cols <- grep("Zscore", colnames(ratio_list))
    for (sample_index in 1:length(zscore_cols)) {
      zscore_col <- zscore_cols[sample_index]
      # matching intensity column
      int_col <- intensity_cols[sample_index]
      # test on column names
      if (check_same_samplename(colnames(ratio_list)[int_col], colnames(ratio_list)[zscore_col])) {
        # calculate Z-scores
        ratio_list[, zscore_col] <- (ratio_list[, int_col] - ratio_list[, "Mean_controls"]) / ratio_list[, "SD_controls"]
      }
    }

    # Add rows of the ratio hmdb codes to the data of zscores from the pipeline.
    dims_xls_ratios <- rbind(ratio_list, dims_xls_copy)

    # Edit the DIMS output Zscores of all patients in format:
    # HMDB_code patientname1  patientname2
    names(dims_xls_ratios) <- gsub("HMDB.code", "HMDB_code", names(dims_xls_ratios))
    names(dims_xls_ratios) <- gsub("HMDB.name", "HMDB_name", names(dims_xls_ratios))

    # for debugging:
    write.table(dims_xls_ratios, file = paste0(output_dir, "/ratios.txt"), sep = "\t")

    # Select only the cols with zscores of the patients
    zscore_patients <- dims_xls_ratios[, c(1, 2, zscore_cols[grep("P", colnames(dims_xls_ratios)[zscore_cols])])]
    # Select only the cols with zscores of the controls
    zscore_controls <- dims_xls_ratios[, c(1, 2, zscore_cols[grep("C", colnames(dims_xls_ratios)[zscore_cols])])]

  }

  #### STEP 4: Run the IEM algorithm #########
  # in: algorithm, file_expected_biomarkers_iem, zscore_patients ||| out: prob_score (+file)
  # algorithm taken from DOI: 10.3390/ijms21030979

  if (algorithm == 1) {
    # Load data
    cat(paste0("\nloading expected file:\n ->  ", file_expected_biomarkers_iem, "\n"))
    expected_biomarkers <- read.csv(file_expected_biomarkers_iem, sep = ";", stringsAsFactors = FALSE)
    # modify column names
    names(expected_biomarkers) <- gsub("HMDB.code", "HMDB_code",  names(expected_biomarkers))
    names(expected_biomarkers) <- gsub("Metabolite", "HMDB_name", names(expected_biomarkers))

    # prepare dataframe scaffold rank_patients
    rank_patients <- zscore_patients
    # Fill df rank_patients with the ranks for each patient
    for (patient_index in 3:ncol(zscore_patients)) {
      # number of positive zscores in patient
      pos <- sum(zscore_patients[, patient_index] > 0)
      # sort the column on zscore; NB: this sorts the entire object, not just one column
      rank_patients <- rank_patients[order(-rank_patients[patient_index]), ]
      # Rank all positive zscores highest to lowest
      rank_patients[1:pos, patient_index] <- as.numeric(ordered(-rank_patients[1:pos, patient_index]))
      # Rank all negative zscores lowest to highest
      rank_patients[(pos + 1):nrow(rank_patients), patient_index] <- as.numeric(ordered(rank_patients[(pos + 1):
                            nrow(rank_patients), patient_index]))
    }

    # Calculate metabolite score, using the dataframes with only values, and later add the cols without values (1&2).
    expected_zscores <- merge(x = expected_biomarkers, y = zscore_patients, by.x = c("HMDB_code"), by.y = c("HMDB_code"))
    expected_zscores_original <- expected_zscores

    # determine which columns contain Z-scores and which contain disease info
    select_zscore_cols <- grep("_Zscore", colnames(expected_zscores))
    select_info_cols <- 1:(min(select_zscore_cols) - 1)
    # set some zscores to zero
    select_incr_indisp <- which(expected_zscores$Change == "Increase" & expected_zscores$Dispensability == "Indispensable")
    expected_zscores[select_incr_indisp, select_zscore_cols] <- lapply(expected_zscores[select_incr_indisp, 
                      select_zscore_cols], function(x) ifelse (x <= 1.6, 0, x))
    select_decr_indisp <- which(expected_zscores$Change == "Decrease" & expected_zscores$Dispensability == "Indispensable")
    expected_zscores[select_decr_indisp, select_zscore_cols] <- lapply(expected_zscores[select_decr_indisp, 
                      select_zscore_cols], function(x) ifelse (x >= -1.2, 0, x))

    # calculate rank score:
    expected_ranks <- merge(x = expected_biomarkers, y = rank_patients, by.x = c("HMDB_code"), by.y = c("HMDB_code"))
    rank_scores <- expected_zscores[order(expected_zscores$HMDB_code), select_zscore_cols] / 
            (expected_ranks[order(expected_ranks$HMDB_code), select_zscore_cols] * 0.9)
    # combine disease info with rank scores
    expected_metabscore <- cbind(expected_ranks[order(expected_zscores$HMDB_code), select_info_cols], rank_scores)

    # multiply weight score and rank score
    weight_score <- expected_zscores
    weight_score[, select_zscore_cols] <- expected_metabscore$Total_Weight * expected_metabscore[, select_zscore_cols]

    # sort table on Disease and Absolute_Weight
    weight_score <- weight_score[order(weight_score$Disease, weight_score$Absolute_Weight, decreasing = TRUE), ]

    # select columns to check duplicates
    dup <- weight_score[, c("Disease", "M.z")]
    uni <- weight_score[!duplicated(dup) | !duplicated(dup, fromLast = FALSE), ]

    # calculate probability score
    prob_score <-  aggregate(uni[, select_zscore_cols], uni["Disease"], sum)

    # list of all diseases that have at least one metabolite Zscore at 0
    for (patient_index in 2:ncol(prob_score)) {
      patient_zscore_colname <- colnames(prob_score)[patient_index]
      matching_colname_expected <- which(colnames(expected_zscores) == patient_zscore_colname)
      # determine which Zscores are 0 for this patient
      zscores_zero <- which(expected_zscores[, matching_colname_expected] == 0)
      # get Disease for these
      disease_zero <- unique(expected_zscores[zscores_zero, "Disease"])
      # set the probability score of these diseases to 0
      prob_score[which(prob_score$Disease %in% disease_zero), patient_index] <- 0
    }

    # determine disease rank per patient
    disease_rank <- prob_score
    # rank diseases in decreasing order
    disease_rank[2:ncol(disease_rank)] <- lapply(2:ncol(disease_rank), function(x) 
                  as.numeric(ordered(-disease_rank[1:nrow(disease_rank), x])))
    # modify column names, Zscores have now been converted to probability scores
    colnames(prob_score) <- gsub("_Zscore", "_prob_score", colnames(prob_score))
    colnames(disease_rank) <- gsub("_Zscore", "", colnames(disease_rank))

    # Create conditional formatting for output Excel sheet. Colors according to values.
    wb <- createWorkbook()
    addWorksheet(wb, "Probability Scores")
    writeData(wb, "Probability Scores", prob_score)
    conditionalFormatting(wb, "Probability Scores", cols = 2:ncol(prob_score), rows = 1:nrow(prob_score), 
        type = "colourScale", style = c("white", "#FFFDA2", "red"), rule = c(1, 10, 100))
    saveWorkbook(wb, file = paste0(output_dir, "/dIEM_algoritme_output_", run_name, ".xlsx"), overwrite = TRUE)
    # check whether prob_score df exists and has expected dimensions.
    if (exists("expected_biomarkers") & (length(disease_rank) == length(prob_score))) {
      cat("\n### Step 4 # Running the IEM algorithm is done.\n\n")
    } else {
      cat("\n**** Error: Could not run IEM algorithm. Check if path to expected_biomarkers csv-file is correct. \n")
    }

    rm(wb)
  }

  #### STEP 5: Make violin plots #####
  # in: algorithm / zscore_patients, violin, nr_contr, nr_pat, Data, path_textfiles, zscore_cutoff, xaxis_cutoff, 
  # top_diseases, top_metab, output_dir ||| out: pdf file, Helix csv file

  if (violin == 1) {

    # preparation
    zscore_patients_copy <- zscore_patients
    # for robust scaler, rename Z-score columns
    colnames(zscore_patients) <- gsub("_RobustZscore", "_Zscore", colnames(zscore_patients))
    colnames(zscore_controls) <- gsub("_RobustZscore", "_Zscore", colnames(zscore_controls))
    colnames(zscore_patients) <- gsub("_Zscore", "", colnames(zscore_patients))
    colnames(zscore_controls) <- gsub("_Zscore", "", colnames(zscore_controls))

    # Make patient list for violin plots
    patient_list <- names(zscore_patients)[-c(1, 2)]

    # from table expected_biomarkers, choose selected columns
    select_columns <- c("Disease", "HMDB_code", "HMDB_name")
    #select_col_nrs <- which(colnames(expected_biomarkers) %in% select_columns)
    expected_biomarkers_select <- expected_biomarkers %>% select(all_of(select_columns))    
    # remove duplicates
    expected_biomarkers_select <- expected_biomarkers_select[!duplicated(expected_biomarkers_select[, c(1, 2)]), ]

    # load file with explanatory information to be included in PDF.
    explanation <- readLines(file_explanation)

    # first step: normal violin plots
    # Find all text files in the given folder, which contain metabolite lists of which
    # each file will be a page in the pdf with violin plots.
    # Make a PDF file for each of the categories in metabolite_dirs
    metabolite_dirs <- list.files(path = path_metabolite_groups, full.names = FALSE, recursive = FALSE)
    for (metabolite_dir in metabolite_dirs) {
      # create a directory for the output PDFs
      pdf_dir <- paste(output_dir, metabolite_dir, sep = "/")
      dir.create(pdf_dir, showWarnings = FALSE)
      cat("making plots in category:", metabolite_dir, "\n")

      # get a list of all metabolite files
      metabolite_files <- list.files(path = paste(path_metabolite_groups, metabolite_dir, sep = "/"), 
            pattern = "*.txt", full.names = FALSE, recursive = FALSE)
      # put all metabolites into one list
      metab_list_all <- list()
      metab_list_names <- c()
      cat("making plots from the input files:")
      # open the text files and add each to a list of dataframes (metab_list_all)
      for (file_index in seq_along(metabolite_files)) {
        infile <- metabolite_files[file_index]
        metab_list <- read.table(paste(path_metabolite_groups, metabolite_dir, infile, sep = "/"), 
              sep = "\t", header = TRUE, quote = "")
        # put into list of all lists
        metab_list_all[[file_index]] <- metab_list
        metab_list_names <- c(metab_list_names, strsplit(infile, ".txt")[[1]][1])
        cat(paste0("\n", infile))
      }
      # include list of classes in metabolite list
      names(metab_list_all) <- metab_list_names

      # prepare list of metabolites; max nr_plots_perpage on one page
      metab_interest_sorted <- prepare_data(metab_list_all, zscore_patients)
      metab_interest_controls <- prepare_data(metab_list_all, zscore_controls)
      metab_perpage <- prepare_data_perpage(metab_interest_sorted, metab_interest_controls, nr_plots_perpage, nr_pat, nr_contr)

      # for Diagnostics metabolites to be saved in Helix
      if(grepl("Diagnost", pdf_dir)) {
        # get table that combines DIMS results with stofgroepen/Helix table
        dims_helix_table <- get_patient_data_to_helix(metab_interest_sorted, metab_list_all)
        
        # check if run contains Diagnostics patients (e.g. "P2024M"), not for research runs
        if(any(is_diagnostic_patient(dims_helix_table$Patient))){
          # get output file for Helix
          output_helix <- output_for_helix(protocol_name, dims_helix_table)
          # write output to file
          path_helixfile <- paste0(output_dir, "/output_Helix_", run_name,".csv")
          write.csv(output_helix, path_helixfile, quote = F, row.names = F)
        }
      }
      
      # make violin plots per patient
      for (pt_nr in 1:length(patient_list)) {
        pt_name <- patient_list[pt_nr]
        # for category Diagnostics, make list of metabolites that exceed alarm values for this patient
        # for category Other, make list of top highest and lowest Z-scores for this patient
        if (grepl("Diagnost", pdf_dir)) {
          top_metab_pt <- prepare_alarmvalues(pt_name, dims_helix_table)
        } else {
          top_metab_pt <- prepare_toplist(pt_name, zscore_patients)
        }

        # generate normal violin plots
        create_violin_plots(pdf_dir, pt_name, metab_perpage, top_metab_pt)

      }

    }

    # Second step: dIEM plots in separate directory
    diem_plot_dir <- paste(output_dir, "dIEM_plots", sep = "/")
    dir.create(diem_plot_dir)

    # Select the metabolites that are associated with the top highest scoring IEM, for each patient
    # disease_rank is from step 4: the dIEM algorithm. The lower the value, the more likely.
    for (pt_nr in 1:length(patient_list)) {
      pt_name <- patient_list[pt_nr]
      # get top diseases for this patient
      pt_colnr <- which(colnames(disease_rank) == pt_name)
      pt_top_indices <- which(disease_rank[, pt_colnr] <= top_nr_iem)
      pt_iems <- disease_rank[pt_top_indices, "Disease"]
      pt_top_iems <- pt_prob_score_top_iems <- c()
      for (single_iem in pt_iems) {
        # get the probability score
        prob_score_iem <- prob_score[which(prob_score$Disease == single_iem), pt_colnr]
        # use only diseases for which probability score is above threshold
        if (prob_score_iem >= threshold_iem) {
          pt_top_iems <- c(pt_top_iems, single_iem)
          pt_prob_score_top_iems <- c(pt_prob_score_top_iems, prob_score_iem)
        }
      }

      # prepare data for plotting dIEM violin plots
      # If prob_score_top_iem is an empty list, don't make a plot
      if (length(pt_top_iems) > 0) {
        # Sorting from high to low, both prob_score_top_iems and pt_top_iems.
        pt_prob_score_order <- order(-pt_prob_score_top_iems)
        pt_prob_score_top_iems <- round(pt_prob_score_top_iems, 1)
        pt_prob_score_top_iem_sorted <- pt_prob_score_top_iems[pt_prob_score_order]
        pt_top_iem_sorted <- pt_top_iems[pt_prob_score_order]
        # getting metabolites for each top_iem disease exactly like in metab_list_all
        metab_iem_all <- list()
        metab_iem_names <- c()
        for (single_iem_index in 1:length(pt_top_iem_sorted)) {
          single_iem <- pt_top_iem_sorted[single_iem_index]
          single_prob_score <- pt_prob_score_top_iem_sorted[single_iem_index]
          select_rows <- which(expected_biomarkers_select$Disease == single_iem)
          metab_list <- expected_biomarkers_select[select_rows, ]
          metab_iem_names <- c(metab_iem_names, paste0(single_iem, ", probability score ", single_prob_score))
          metab_list <- metab_list[, -1]
          metab_iem_all[[single_iem_index]] <- metab_list
        }
        # put all metabolites into one list
        names(metab_iem_all) <- metab_iem_names

        # get Zscore information from zscore_patients_copy, similar to normal violin plots
        metab_iem_sorted <- prepare_data(metab_iem_all, zscore_patients_copy)
        metab_iem_controls <- prepare_data(metab_iem_all, zscore_controls)
        # make sure every page has 20 metabolites
        diem_metab_perpage <- prepare_data_perpage(metab_iem_sorted, metab_iem_controls, nr_plots_perpage, nr_pat)
	# add table of metabolites with increased or decreased Z-scores
        top_metab_pt <- prepare_toplist(pt_name, zscore_patients)

        # generate dIEM violin plots
        create_violin_plots(diem_plot_dir, pt_name, diem_metab_perpage, top_metab_pt)

      } else {
        cat(paste0("\n\n**** This patient had no prob_scores higher than ", threshold_iem, ".
                    Therefore, this pdf was not made:\t ", pt_name, "_iem \n"))
      }

    }

  }
}
