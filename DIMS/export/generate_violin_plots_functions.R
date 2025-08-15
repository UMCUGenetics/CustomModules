#' Getting the intensities for calculating ratio Z-scores
#'
#' @param ratios_metabs_df: dataframe with HMDB codes for the ratios (dataframe)
#' @param row_index: index of the row in the ratios_metabs_df (integer)
#' @param intensities_zscore_df: dataframe with intensities for each sample (dataframe)
#' @param fraction_side: either numerator or denominator, which side of the fraction (string)
#' @param intensity_cols: names of the columns that contain the intensities (string)
#'
#' @returns fraction_side_intensity: a vector of intensities (vector of integers)
get_intentities_for_ratios <- function(ratios_metabs_df, row_index, intensities_zscore_df, fraction_side, intensity_cols) {
  fraction_side_hmdb_ids <- ratios_metabs_df[row_index, fraction_side]
  if (grepl("plus", fraction_side_hmdb_ids)) {
    fraction_side_hmdb_id_list <- strsplit(fraction_side_hmdb_ids, "plus")[[1]]
    fraction_side_intensity_list <- intensities_zscore_df %>% filter(HMDB_code %in% fraction_side_hmdb_id_list) %>%
      select(any_of(intensity_cols))
    fraction_side_intensity <- apply(fraction_side_intensity_list, 2, sum)
  } else if(fraction_side_hmdb_ids == "one") {
    fraction_side_intensity <- 1
  } else {
    fraction_side_intensity <- intensities_zscore_df %>% filter(HMDB_code == fraction_side_hmdb_ids) %>%
      select(any_of(intensity_cols))
  }
  return(fraction_side_intensity)
}

#' Get the sample IDs for columns that have Z-score and intensities
#'
#' @param colnames_zscore: vector of sample IDs from the dataframe containing Z-scores (vector of strings)
#' @param intensity_cols: vector of sample IDs form the dataframe containing intensities (vector of strings)
#'
#' @returns: vector of sample IDs that are in both input vectors (vector of strings)
get_zscore_columns <- function(colnames_zscore, intensity_cols) {
  sample_intersect <- intersect(paste0(intensity_cols, "_Zscore"), grep("_Zscore", colnames_zscore, value = TRUE))
  return(sample_intersect)
}

#' Get a list with dataframes for all off the metabolite group in a directory
#'
#' @param metab_group_dir: directory containing txt files with metabolites per group (string)
#'
#' @returns: list with dataframes with info on metabolites (list of dataframes)
get_list_metabolites <- function(metab_group_dir) {
  # get a list of all metabolite files
  metabolite_files <- list.files(metab_group_dir, pattern = "*.txt", full.names = FALSE, recursive = FALSE)
  # put all metabolites into one list
  metab_list_all <- lapply(paste(metab_group_dir, metabolite_files, sep = "/"), 
                           read.table, sep = "\t", header = TRUE, quote = "")
  names(metab_list_all) <- gsub(".txt", "", metabolite_files)
  
  return(metab_list_all)
}

#' Combine patient Z-scores with metabolite info
#'
#' @param metab_list_all: list of dataframes with metabolite information for different stofgroepen (list)
#' @param zscore_df: dataframe with metabolite Z-scores for all patient
#'
#' @return: list of dataframes for each stofgroep with data for each metabolite and patient/control per row
combine_metab_info_zscores <- function(metab_list_all, zscore_df) {
  # remove HMDB_name column and "_Zscore" from column (patient) names
  zscore_df <- zscore_df %>% select(-HMDB_name) %>%
    rename_with(~ str_remove(.x, "_Zscore"), .cols = contains("_Zscore"))
  
  # put data into pages, max 20 violin plots per page in PDF
  metab_interest_sorted <- list()
  
  for (metab_class in names(metab_list_all)) {
    metab_df <- metab_list_all[[metab_class]]
    # Select HMDB_code and HMDB_name columns
    metab_df <- metab_df %>% select(HMDB_code, HMDB_name)
    
    # Change the HMDB_name column so all names have 45 characters
    metab_df <- metab_df %>% mutate(HMDB_name = case_when(
      str_length(HMDB_name) > 45 ~ str_c(str_sub(HMDB_name, 1, 42), "..."),
      str_length(HMDB_name) < 45 ~ str_pad(HMDB_name, 45, side = "right", pad = " "),
      TRUE ~ HMDB_name
    ))
    
    # Join metabolite info with the Z-score dataframe
    metab_interest <- metab_df %>% inner_join(zscore_df, by = "HMDB_code") %>% select(-HMDB_code)
    
    # put the data frame in long format
    metab_interest_melt <- reshape2::melt(metab_interest, id.vars = "HMDB_name", variable.name = "Sample", 
                                          value.name = "Z_score")
    # Add the dataframe sorted on HMDB_name to a list
    metab_interest_sorted[[metab_class]] <- metab_interest_melt
  }
  
  return(metab_interest_sorted)
}

#' Combine patient and control data for each page of the violinplot pdf
#'
#' @param metab_interest_sorted: list of dataframes with data for each metabolite and patient (list)
#' @param metab_interest_contr:  list of dataframes with data for each metabolite and control (list)
#' @param nr_plots_perpage: number of plots per page in the violinplot pdf (integer)
#' @param nr_pat: number of patients (integer)
#' @param nr_contr: number of controls (integer)
#'
#' @return: list of dataframes with metabolite Z-scores for each patient and control,
#'          the length of list is the number of pages for the violinplot pdf (list)
prepare_data_perpage <- function(metab_interest_sorted, metab_interest_contr, nr_plots_perpage, nr_pat, nr_contr) {
  metab_perpage <- list()
  metab_category <- c()
  
  for (metab_class in names(metab_interest_sorted)) {
    # Get the data for patients and controls for the metab_interest_sorted list
    metab_sort_patients_df <- metab_interest_sorted[[metab_class]]
    metab_sort_controls_df <- metab_interest_contr[[metab_class]]
    
    # Calculate the number of pages
    nr_pages <- ceiling(length(unique(metab_sort_patients_df$HMDB_name)) / nr_plots_perpage)
    
    # Get all metabolites and create list with HMDB naames of max nr_plots_perpage long
    metabolites <- unique(metab_sort_patients_df$HMDB_name)
    metabolites_in_chunks <- split(metabolites, ceiling(seq_along(metabolites) / nr_plots_perpage))
    nr_chunks <- length(metabolites_in_chunks)
    
    current_perpage <- lapply(metabolites_in_chunks, function(metab_name) {
      patients_df <- metab_sort_patients_df %>% filter(HMDB_name %in% metab_name)
      controls_df <- metab_sort_controls_df %>% filter(HMDB_name %in% metab_name)
      
      # Combine both dataframes
      combined_df <- rbind(patients_df, controls_df)
      
      # Add empty dummy's to extend the number of metabs to the nr_plots_perpage
      n_missing <- nr_plots_perpage - length(metab_name)
      if (n_missing > 0) {
        dummy_names <- paste0(" ", strrep(" ", seq_len(n_missing)))
        metab_order <- c(metab_name, dummy_names)
      } else {
        metab_order <- metab_name
      }
      attr(combined_df, "y_order") <- rev(metab_order)
      
      return(combined_df)
    })
    # Add new items to main list
    metab_perpage <- append(metab_perpage, current_perpage)
    # create list of page headers
    metab_category <- c(metab_category, paste(metab_class, seq(nr_chunks), sep = "_"))
  }
  # add page headers to list
  names(metab_perpage) <- metab_category
  
  return(metab_perpage)
}

#' Get patient data to be uploaded to Helix
#'
#' @param metab_interest_sorted: list of dataframes with metabolite Z-scores for each sample/patient (list)
#' @param metab_list_all: list of tables with metabolites for Helix and violin plots (list)
#'
#' @return: dataframe with patient data with only metabolites for Helix and violin plots 
#'          with Helix name, high/low Z-score cutoffs
get_patient_data_to_helix <- function(metab_interest_sorted, metab_list_all) {
  # Combine Z-scores of metab groups together
  df_all_metabs_zscores <- bind_rows(metab_interest_sorted)
  
  # Change the Sample column to characters, trim HMDB_name and split HMDB_name in new column
  df_all_metabs_zscores <- df_all_metabs_zscores %>%
    mutate(Sample = as.character(Sample),
           HMDB_name = str_trim(HMDB_name, "right"),
           HMDB_name_split = str_split_fixed(HMDB_name, "nitine;", 2)[, 1])
  
  # Combine stofgroepen
  dims_helix_table <- bind_rows(metab_list_all)
  
  # Filter for Helix metabolites and split HMDB_name column for matching with df_all_metabs_zscores
  dims_helix_table <- dims_helix_table %>%
    filter(Helix == "ja") %>%
    mutate(HMDB_name_split = str_split_fixed(HMDB_name, "nitine;", 2)[, 1]) %>%
    select(HMDB_name_split, Helix_naam, high_zscore, low_zscore)
  
  # Filter DIMS results for metabolites for Helix and combine Helix info
  df_metabs_helix <- df_all_metabs_zscores %>%
    filter(HMDB_name_split %in% dims_helix_table$HMDB_name_split) %>%
    left_join(dims_helix_table, by = join_by(HMDB_name_split)) %>%
    select(HMDB_name, Sample, Z_score, Helix_naam, high_zscore, low_zscore)
  
  return(df_metabs_helix)
}

#' Check for Diagnostics patients with correct patient number (e.g. starting with "P2024M")
#'
#' @param patient_column: a column from dataframe with IDs (character vector)
#'
#' @return: a logical vector with TRUE or FALSE for each element (vector)
is_diagnostic_patient <- function(patient_column) {
  diagnostic_patients <- grepl("^P[0-9]{4}M", patient_column)
  
  return(diagnostic_patients)
}

#' Get the output dataframe for Helix
#'
#' @param protocol_name: protocol name (string)
#' @param df_metabs_helix: dataframe with metabolite Z-scores for patients (dataframe)
#'
#' @return: dataframe with patient metabolite Z-scores in correct format for Helix
output_for_helix <- function(protocol_name, df_metabs_helix) {
  # Remove positive controls
  df_metabs_helix <- df_metabs_helix %>% filter(is_diagnostic_patient(Sample))
  
  # Add 'Vial' column, each patient has unique ID
  df_metabs_helix <- df_metabs_helix %>%
    group_by(Sample) %>%
    mutate(Vial = cur_group_id()) %>%
    ungroup()
  
  # Split patient number into labnummer and Onderzoeksnummer
  df_metabs_helix <- add_lab_id_and_onderzoeksnummer(df_metabs_helix)
  
  # Add column with protocol name
  df_metabs_helix$Protocol <- protocol_name
  
  # Change name Z_score and Helix_naam columns to Amount and Name
  change_columns <- c(Amount = "Z_score", Name = "Helix_naam")
  df_metabs_helix <- df_metabs_helix %>% rename(all_of(change_columns))
  
  # Select only necessary columns and set them in correct order
  df_metabs_helix <- df_metabs_helix %>%
    select(c(Vial, labnummer, Onderzoeksnummer, Protocol, Name, Amount))
  
  # Remove duplicate patient-metabolite combinations ("leucine + isoleucine + allo-isoleucin_Z-score" is added 3 times)
  df_metabs_helix <- df_metabs_helix %>%
    group_by(Onderzoeksnummer, Name) %>%
    distinct() %>%
    ungroup()
  
  return(df_metabs_helix)
}

#' Adding labnummer and Onderzoeksnummer to a dataframe
#'
#' @param df_metabs_helix: dataframe with patient data to be uploaded to Helix
#'
#' @return: dataframe with added labnummer and Onderzoeksnummer columns
add_lab_id_and_onderzoeksnummer <- function(df_metabs_helix) {
  # Split patient number into labnummer and Onderzoeksnummer
  for (row in 1:nrow(df_metabs_helix)) {
    df_metabs_helix[row, "labnummer"] <- gsub("^P|\\.[0-9]*", "", df_metabs_helix[row, "Sample"])
    labnummer_split <- strsplit(as.character(df_metabs_helix[row, "labnummer"]), "M")[[1]]
    df_metabs_helix[row, "Onderzoeksnummer"] <- paste0("MB", labnummer_split[1], "/", labnummer_split[2])
  }
  
  return(df_metabs_helix)
}

#' Create a dataframe with all metabolites that exceed the min and max Z-score cutoffs
#'
#' @param patient_name: patient code (string)
#' @param dims_helix_table: dataframe with metabolite Z-scores for each patient and Helix info (dataframe)
#'
#' @return: dataframe with metabolites that exceed the min and max Z-score cutoffs for the selected patient
prepare_alarmvalues <- function(patient_name, dims_helix_table) {
  # extract data for patient of interest (patient_name)
  patient_metabs_helix <- dims_helix_table %>%
    filter(Sample == patient_name) %>%
    mutate(Z_score = round(Z_score, 2))
  
  patient_high_df <- patient_metabs_helix %>% filter(Z_score > high_zscore)
  patient_low_df <- patient_metabs_helix %>% filter(Z_score < low_zscore)
  
  if (nrow(patient_high_df) > 0 | nrow(patient_low_df) > 0) {
    # sort tables on zscore
    patient_high_df <- patient_high_df %>% arrange(desc(Z_score)) %>% select(c(HMDB_name, Z_score))
    patient_low_df <- patient_low_df %>% arrange(Z_score) %>% select(c(HMDB_name, Z_score))
  }
  # add lines for increased, decreased
  extra_line1 <- c("Increased", "")
  extra_line2 <- c("Decreased", "")
  
  # combine the two lists
  top_metab_patient <- rbind(extra_line1, patient_high_df, extra_line2, patient_low_df)
  
  # remove row names
  rownames(top_metab_patient) <- NULL
  # change column names for display
  colnames(top_metab_patient) <- c("Metabolite", "Z-score")
  
  return(top_metab_patient)
}

#' Create a dataframe with the top 20 highest and top 10 lowest metabolites per patient
#'
#' @param pt_name: patient code (string)
#' @param zscore_patients: dataframe with metabolite Z-scores per patient (dataframe)
#' @param top_highest: the number of metabolites with the highest Z-score to display in the table (numeric) 
#' @param top_lowest: the number of metabolites with the lowest Z-score to display in the table (numeric)
#'
#' @return: dataframe with 30 metabolites and Z-scores (dataframe)
prepare_toplist <- function(patient_id, zscore_patients) {
  top_highest <- 20
  top_lowest <- 10
  patient_df <- zscore_patients %>%
    select(HMDB_code, HMDB_name, !!sym(patient_id)) %>%
    arrange(!!sym(patient_id))
  
  # Get lowest Zscores
  patient_df_low <- patient_df[1:top_lowest, ]
  patient_df_low <- patient_df_low %>% mutate(across(!!sym(patient_id), ~ round(.x ,2)))
  
  # Get highest Zscores
  patient_df_high <- patient_df[nrow(patient_df):(nrow(patient_df) - top_highest + 1), ]
  patient_df_high <- patient_df_high %>% mutate(across(!!sym(patient_id), ~ round(.x ,2)))
  
  # add lines for increased, decreased
  extra_line1 <- c("Increased", "", "")
  extra_line2 <- c("Decreased", "", "")
  top_metab_pt <- rbind(extra_line1, patient_df_high, extra_line2, patient_df_low)
  # remove row names
  rownames(top_metab_pt) <- NULL
  
  # change column names for display
  colnames(top_metab_pt) <- c("HMDB_ID", "Metabolite", "Z-score")
  
  return(top_metab_pt)
}

#' Create a pdf with table with metabolites and violin plots
#'
#' @param pdf_dir: location where to save the pdf file (string)
#' @param patient_id: patient id (string)
#' @param metab_perpage: list of dataframes, each dataframe contains data for a page in de pdf (list)
#' @param top_metab_pt: dataframe with increased and decreased metabolites for this patient (dataframe)
#' @param explanation: text that explains the violin plots and the pipeline version (string)
create_pdf_violin_plots <- function(pdf_dir, patient_id, metab_perpage, top_metab_pt, explanation) {
  # set parameters for plots
  plot_height <- 9.6
  plot_width <- 6
  
  # patient plots, create the PDF device
  patient_id_sub <- patient_id
  suffix <- ""
  if (grepl("Diagnostics", pdf_dir) & is_diagnostic_patient(patient_id)) {
    prefix <- "MB"
    suffix <- "_DIMS_PL_DIAG"
    # substitute P and M in P2020M00001 into right format for Helix
    patient_id_sub <- gsub("[PM]", "", patient_id)
    patient_id_sub <- gsub("\\..*", "", patient_id_sub)
  } else if (grepl("Diagnostics", pdf_dir)) {
    prefix <- "Dx_"
  } else if (grepl("IEM", pdf_dir)) {
    prefix <- "IEM_"
  } else {
    prefix <- "R_"
  }
  
  pdf(paste0(pdf_dir, "/", prefix, patient_id_sub, suffix, ".pdf"),
      onefile = TRUE,
      width = plot_width,
      height = plot_height)
  
  # page headers:
  page_headers <- names(metab_perpage)
  
  # put table into PDF file, if not empty
  if (!is.null(dim(top_metab_pt))) {
    max_rows_per_page <- 35
    total_rows <- nrow(top_metab_pt)
    number_of_pages <- ceiling(total_rows / max_rows_per_page)
    
    # get the names and numbers in the table aligned
    table_theme <- ttheme_default(core = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 6)),
                                  colhead = list(fg_params = list(fontsize = 8, fontface = "bold")))
    
    for (page in seq(number_of_pages)) {
      start_row <- (page - 1) * max_rows_per_page + 1
      end_row <- min(page * max_rows_per_page, total_rows)
      page_data <- top_metab_pt[start_row:end_row, ]
      
      table_grob <- tableGrob(page_data, theme = table_theme, rows = NULL)
      
      grid.arrange(
        table_grob,
        top = paste0("Top deviating metabolites for patient: ", patient_id)
      )
    }
  }
  
  # violin plots
  for (metab_class in names(metab_perpage)) {
    # extract list of metabolites to plot on a page
    metab_zscores_df <- metab_perpage[[metab_class]]
    # extract original data for patient of interest (pt_name) before cut-offs
    patient_zscore_df <- metab_zscores_df %>% filter(Sample == patient_id)
    
    # Remove patient column and change Z-score. If under -5 to -5 and if above 20 to 20.    
    metab_zscores_df <- metab_zscores_df %>%
      filter(Sample != patient_id) %>%
      mutate(Z_score = pmin(pmax(Z_score, -5), 20))
    
    # subtitle per page
    sub_perpage <- gsub("_", " ", metab_class)
    # for IEM plots, put subtitle on two lines
    sub_perpage <- gsub("probability", "\nprobability", sub_perpage)
    
    # draw violin plot. 
    ggplot_object <- create_violin_plot(metab_zscores_df, patient_zscore_df, sub_perpage, patient_id)
    
    suppressWarnings(print(ggplot_object))
  }
  
  # add explanation of violin plots, version number etc.
  plot(NA, xlim = c(0, 5), ylim = c(0, 5), bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  if (length(explanation) > 0) {
    text(0.2, 5, explanation[1], pos = 4, cex = 0.8)
    for (line_index in 2:length(explanation)) {
      text_y_position <- 5 - (line_index * 0.2)
      text(-0.2, text_y_position, explanation[line_index], pos = 4, cex = 0.5)
    }
  }
  
  # close the PDF file
  dev.off()
}

#' Create violin plots
#'
#' @param metab_zscores_df: dataframe with Z-scores for all samples (dataframe)
#' @param patient_zscore_df: dataframe with Z-scores for the specified patient (dataframe)
#' @param sub_perpage: subtitle of the page (string) 
#' @param patient_id: the patient id of the selected patient (string)
#'
#' @returns ggpplot_object: a violin plot of metabolites that highlights the selected patient (ggplot object)
create_violin_plot <- function(metab_zscores_df, patient_zscore_df, sub_perpage, patient_id) {
  fontsize <- 1
  circlesize <- 0.8
  # Set colors for the violinplot: green, blue, blue/purple, purple, orange, red
  colors_plot <- c("#22E4AC", "#00B0F0", "#504FFF", "#A704FD", "#F36265", "#DA0641")
  
  y_order <- attr(metab_zscores_df, "y_order")
  metab_zscores_df$HMDB_name <- rev(factor(metab_zscores_df$HMDB_name, levels = rev(y_order)))
  patient_zscore_df$HMDB_name <- rev(factor(patient_zscore_df$HMDB_name, levels = rev(y_order)))
  
  ggplot_object <- ggplot(metab_zscores_df, aes(x = Z_score, y = HMDB_name)) +
    # Make violin plots
    geom_violin(scale = "width", na.rm = TRUE) +
    # Add Z-score for the selected patient, shape=22 gives square for patient of interest
    geom_point(data = patient_zscore_df, aes(color = Z_score),
               size = 3.5 * circlesize, shape = 22, fill = "white", na.rm = TRUE) +
    # Add the Z-score at the right side of the plot
    geom_text(
      data = patient_zscore_df,
      aes(16, label = paste0("Z=", round(Z_score, 2))),
      hjust = "left", vjust = +0.2, size = 3, na.rm = TRUE) +
    # Set colour for the Z-score of the selected patient
    scale_fill_gradientn(
      colors = colors_plot, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar",
      aesthetics = "colour"
    ) +
    # Add labels to the axis
    labs(x = "Z-scores", y = "Metabolites", subtitle = sub_perpage, color = "z-score") +
    # Add a title to the page
    ggtitle(label = paste0("Results for patient ", patient_id)) +
    # Set theme: size and font type of y-axis labels, remove legend and make the 
    theme(
      axis.text.y = element_text(family = "Courier", size = 6),
      legend.position = "none",
      plot.caption = element_text(size = rel(fontsize))
    ) +
    # Set y-axis to set order
    scale_y_discrete(limits = y_order) +
    # Limit the x-axis to between -5 and 20
    xlim(-5, 20) +
    # Set grey vertical lines at -2 and 2
    geom_vline(xintercept = c(-2, 2), col = "grey", lwd = 0.5, lty = 2)
  
  
  return(ggplot_object)
}

#' Run the dIEM algorithm (DOI: 10.3390/ijms21030979)
#'
#' @param expected_biomarkers_df: table with information for HMDB codes about IEMs (dataframe) 
#' @param zscore_patients: dataframe containing Z-scores for patient (dataframe)
#'
#' @returns probability_score: a dataframe with probability scores for IEMs for each patient (dataframe)
run_diem_algorithm <- function(expected_biomarkers_df, zscore_patients_df, sample_cols) {
  # Rank the metabolites for each patient individually
  ranking_patients <- zscore_patients_df %>% 
    mutate(across(-c(HMDB_code, HMDB_name), rank_patient_zscores))
  
  ranking_patients <- merge(x = expected_biomarkers_df, y = ranking_patients,
                            by.x = c("HMDB_code"), by.y = c("HMDB_code"))
  
  zscore_expected_df <- merge(x = expected_biomarkers_df, y = zscore_patients_df,
                              by.x = c("HMDB_code"), by.y = c("HMDB_code"))
  
  # Change Z-score to zero for specific cases
  zscore_expected_df <- zscore_expected_df %>% mutate(across(
    all_of(sample_cols),
    ~ case_when(
      Change == "Increase" & Dispensability == "Indispensable" & .x <= 1.6 ~ 0,
      Change == "Decrease" & Dispensability == "Indispensable" & .x >= -1.2 ~ 0,
      TRUE ~ .x
    )
  ))
  
  # Sort both dataframes on HMDB_code for calculating the metabolite score
  zscore_expected_df <- zscore_expected_df[order(zscore_expected_df$HMDB_code), ]
  ranking_patients <- ranking_patients[order(ranking_patients$HMDB_code), ]
  
  # Set up dataframe for the metabolite score, copy zscore_expected_df for biomarker info
  metabolite_score_info <- zscore_expected_df
  # Calculate metabolite score: Z-score/(Rank * 0.9)
  metabolite_score_info[sample_cols] <- zscore_expected_df[sample_cols] / (ranking_patients[sample_cols] * 0.9)
  
  # Calculate the weighted score: metabolite_score * Total_Weight
  metabolite_weight_score <- metabolite_score_info %>% 
    mutate(across(
      all_of(sample_cols),
      ~ .x * Total_Weight
    )) %>%
    arrange(desc(Disease), desc(Absolute_Weight))
  
  # Calculate the probability score for each disease - Mz combination
  probability_score <- metabolite_weight_score %>%
    filter(
      !duplicated(select(., Disease, M.z)) |
      !duplicated(select(., Disease, M.z), fromLast = FALSE)
    ) %>%
    group_by(Disease) %>% 
    summarise(across(all_of(sample_cols), sum), .groups = "drop")
  
  
  # Set probability score to 0 for Z-scores == 0
  for (sample_col in sample_cols) {
    # Get indexes of Zscore that equal 0
    zscores_zero_idx <- which(zscore_expected_df[[sample_col]] == 0)
    # Get diseases that have a Zscore of 0
    diseases_zero <- unique(zscore_expected_df[zscores_zero_idx, "Disease"])
    # Set probabilty of these diseases to 0
    probability_score[probability_score$Disease %in% diseases_zero, sample_col] <- 0
  }
  
  colnames(probability_score) <- gsub("_Zscore", "_prob_score", colnames(probability_score))
  
  return(probability_score)
}

#' Ranking Z-scores for a patient, separate for positive and negative Z-scores
#'
#' @param zscore_col: vector with Z-scores for a single patient (vector of integers)
#'
#' @returns ranking: a vector of the ranking of the Z-scores (vector of integers)
rank_patient_zscores <- function(zscore_col) {
  # Create ranking column with default NA values
  ranking <- rep(NA_real_, length(zscore_col))
  
  # Get indexes for negative and positive rows
  neg_indexes <- which(zscore_col <= 0)
  pos_indexes <- which(zscore_col > 0)
  
  # Rank the negative and positive Zscores
  ranking[neg_indexes] <- dense_rank(zscore_col[neg_indexes])
  ranking[pos_indexes] <- dense_rank(-zscore_col[pos_indexes])
  
  return(ranking)
}

#' Save the probability score dataframe as an Excel file
#'
#' @param probability_score: a dataframe containing probability scores for each patient (dataframe)
#' @param output_dir: location where to save the Excel file (string)
#' @param run_name: name of the run, for the file name (string)
save_prob_scores_to_Excel <- function(probability_score, output_dir, run_name) {
  # Create conditional formatting for output Excel sheet. Colors according to values.
  wb <- createWorkbook()
  addWorksheet(wb, "Probability Scores")
  writeData(wb, "Probability Scores", probability_score)
  conditionalFormatting(wb, "Probability Scores", cols = 2:ncol(probability_score), rows = 1:nrow(probability_score), 
                        type = "colourScale", style = c("white", "#FFFDA2", "red"), rule = c(1, 10, 100))
  saveWorkbook(wb, file = paste0(output_dir, "/dIEM_algoritme_output_", run_name, ".xlsx"), overwrite = TRUE)
  rm(wb)
}
