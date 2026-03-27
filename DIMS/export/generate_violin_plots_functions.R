#' Preparing the intensities and Z-score dataframe.
#' Certain columns are removed, the HMDB_code and HMDB_name column are moved forward,
#' the avg_ctrls and sd_ctrls columns are renamed and the column type of all columns containing numbers
#' is changed to numeric.
#'
#' @param intensities_zscore_df: dataframe with intensities, Z-scores and metabolite information for all samples
#'
#' @returns intensities_zscore_df: a dataframe containing intensities, Z-scores, HMDB IDs, HMDB names and
#' the mean and average of all controls
prepare_intensities_zscore_df <- function(intensities_zscore_df) {
  intensities_zscore_df <- intensities_zscore_df %>%
    select(-c(
      plots, HMDB_name_all, HMDB_ID_all, sec_HMDB_ID, HMDB_key, sec_HMDB_ID_rlvnc, name,
      relevance, descr, origin, fluids, tissue, disease, pathway, nr_ctrls
    )) %>%
    relocate(c(HMDB_code, HMDB_name)) %>%
    rename(mean_controls = avg_ctrls, sd_controls = sd_ctrls) %>%
    mutate(across(!c(HMDB_name, HMDB_code), as.numeric))
  return(intensities_zscore_df)
}

#' Get all column names containing a specific prefix.
#'
#' @param dataframe: dataframe containing multiple columns with Z-scores
#' @param prefix: a string of a prefix to be searched in the column names, e.g. "P" or "C".
#'
#' @returns sample_colnames: a vector of column names all containing the prefix.
get_colnames_by_prefix <- function(dataframe, prefix) {
  sample_colnames <- grep(paste0("^", prefix), colnames(dataframe), value = TRUE)
  return(sample_colnames)
}

#' Remove the suffix from a vector of names
#'
#' @param vector_names: vector containing names with or without a suffix
#' @param suffix: string containing the suffix to be removed
#'
#' @returns names_no_suffix: a vector of unique names without the suffix
remove_suffix_from_items <- function(vector_names, suffix) {
  names_no_suffix <- unique(gsub("_Zscore", "", vector_names))
  return(names_no_suffix)
}

#' Add Zscores for multiple ratios to the dataframe
#'
#' @param outlist: dataframe containing intensities and Z-scores for all controls and patients
#' @param metabolites_ratios_df: dataframe containing numerators and denominators for all ratios
#' @param all_sample_ids: vector of sample IDS, controls and patients
#'
#' @returns intensities_zscore_ratios_df: dataframe containing intensities and Z-scores for all controls and patients
#' for all metabolites and ratios
add_zscores_ratios_to_df <- function(outlist, metabolites_ratios_df, all_sample_ids) {
  intensities_zscores_df <- prepare_intensities_zscore_df(outlist)

  # calculate Z-scores for the ratios
  zscore_ratios_df <- calculate_zscore_ratios(metabolites_ratios_df, intensities_zscores_df, all_sample_ids)
  intensities_zscore_ratios_df <- rbind(intensities_zscores_df, zscore_ratios_df)

  return(intensities_zscore_ratios_df)
}

#' Calculate Z-scores for ratios
#'
#' @param metabolites_ratios_df: dataframe containing numerators and denominators for all ratios
#' @param intensities_zscores_df: dataframe containing intensities and Z-scores for all controls and patients
#' @param intensity_col_names: vector of sample IDS, controls and patients
#'
#' @returns zscore_ratios_df: dataframe containing Z-scores for all ratios for all samples
calculate_zscore_ratios <- function(metabolites_ratios_df, intensities_zscores_df, intensity_col_names) {
  zscore_ratios_df <- data.frame(matrix(
    ncol = ncol(intensities_zscores_df),
    nrow = nrow(metabolites_ratios_df)
  ))
  colnames(zscore_ratios_df) <- colnames(intensities_zscores_df)

  # put HMDB info into first two columns of ratio_zscore_df
  zscore_ratios_df$HMDB_code <- metabolites_ratios_df$HMDB.code
  zscore_ratios_df$HMDB_name <- metabolites_ratios_df$Ratio_name

  intensity_cols_index <- which(colnames(zscore_ratios_df) %in% intensity_col_names)
  for (row_index in seq_len(nrow(metabolites_ratios_df))) {
    # Get a list of intensities for the numerator
    numerator_intensities <- get_intensities_fraction_side(
      metabolites_ratios_df,
      row_index,
      intensities_zscores_df,
      "HMDB_numerator",
      intensity_col_names
    )
    # Get a list of intensities for the denominator
    denominator_intensities <- get_intensities_fraction_side(
      metabolites_ratios_df,
      row_index,
      intensities_zscores_df,
      "HMDB_denominator",
      intensity_col_names
    )
    # calculate the intensity ratio for each sample
    zscore_ratios_df[row_index, intensity_cols_index] <- log2(numerator_intensities / denominator_intensities)
  }

  control_intensities_cols_index <- grep("^C[^_]*$", colnames(intensities_zscores_df), perl = TRUE)
  # Calculate means and SD's of the calculated ratios for Controls
  zscore_ratios_df[, "mean_controls"] <- apply(zscore_ratios_df[, control_intensities_cols_index], 1, mean)
  zscore_ratios_df[, "sd_controls"] <- apply(zscore_ratios_df[, control_intensities_cols_index], 1, sd)

  # Calculate Zscores for the ratios
  samples_zscore_columns <- get_sample_ids_with_zscores(colnames(intensities_zscores_df), intensity_col_names)
  intensity_ratios_df <- zscore_ratios_df[, intensity_col_names]
  mean_ratios_controls <- zscore_ratios_df[, "mean_controls"]
  sd_ratios_controls <- zscore_ratios_df[, "sd_controls"]

  zscore_ratios_df[, samples_zscore_columns] <- (intensity_ratios_df - mean_ratios_controls) / sd_ratios_controls

  return(zscore_ratios_df)
}

#' Make and save violin plots for each patient in a PDF
#'
#' @param zscore_patients_df: dataframe with Z-scores for all patient samples
#' @param zscore_controls_df: dataframe with Z-scores for all control samples
#' @param path_metabolite_groups: string containing the path for the metabolite groups directories
#' @param nr_plots_perpage: integer containing the number of metabolites on a plot per page
#' @param number_of_samples: list containing the number of patient and control samples
#' @param run_name: string containing the run name
#' @param protocol_name: string containing the protocol name
#' @param explanation_violin_plot: vector of strings containing the explanation of the violin plots
#' @param number_of_metabolites: list containing the number of metabolites for the top and lowest table
make_and_save_violin_plot_pdfs <- function(
    zscore_patients_df,
    zscore_controls_df,
    path_metabolite_groups,
    nr_plots_perpage,
    number_of_samples,
    run_name,
    protocol_name,
    explanation_violin_plot,
    number_of_metabolites) {
  # Get all patient IDs
  patient_col_names <- remove_suffix_from_items(get_colnames_by_prefix(zscore_patients_df, "P"), "_Zscore")
  # get all files from metabolite_groups directory
  metabolite_dirs <- list.files(path = path_metabolite_groups, full.names = FALSE, recursive = FALSE)
  for (metabolite_dir in metabolite_dirs) {
    # create a directory for the output PDFs
    pdf_dir <- paste0("./", metabolite_dir)
    dir.create(pdf_dir, showWarnings = FALSE)

    metab_list_all <- get_list_dataframes_from_dir(paste(path_metabolite_groups, metabolite_dir, sep = "/"))
    metab_interest_patients <- merge_metabolite_info_zscores(metab_list_all, zscore_patients_df)
    metab_interest_controls <- merge_metabolite_info_zscores(metab_list_all, zscore_controls_df)
    metab_perpage <- get_data_per_metabolite_class(
      metab_interest_patients,
      metab_interest_controls,
      nr_plots_perpage,
      number_of_samples$patients,
      number_of_samples$controls
    )

    # for Diagnostics metabolites to be saved in Helix
    if (grepl("Diagnost", pdf_dir)) {
      # get table that combines DIMS results with metabolite classes/Helix table
      dims_helix_table <- prepare_helix_patient_data(metab_interest_patients, metab_list_all)
      # check if run contains diagnostic patients (e.g. "P2024M")
      if (any(is_diagnostic_patients(dims_helix_table$Sample))) {
        # transform dataframe for Helix output
        output_helix <- transform_metab_df_to_helix_df(protocol_name, dims_helix_table)
        # save the DIMS Helix dataframe
        path_helixfile <- paste0("./output_Helix_", run_name, ".csv")
        write.csv(output_helix, path_helixfile, quote = FALSE, row.names = FALSE)
      }
    }

    # make violin plots per patient
    for (patient_id in patient_col_names) {
      if (grepl("Diagnost", pdf_dir)) {
        # make list of metabolites that exceed alarm values for this patient
        top_metabs_patient <- get_top_metabolites_df(patient_id, dims_helix_table)
      } else {
        # make list of top highest and lowest Z-scores for this patient
        top_metabs_patient <- prepare_toplist(
          patient_id,
          zscore_patients_df,
          number_of_metabolites$highest,
          number_of_metabolites$lowest
        )
      }
      # generate normal violin plots
      create_pdf_violin_plots(pdf_dir, patient_id, metab_perpage, top_metabs_patient, explanation_violin_plot)
    }
  }
}

#' Get a list with dataframes for all off the metabolite group in a directory
#'
#' @param dir_with_subdirs: directory containing txt files with metabolites per group (string)
#'
#' @returns list_of_dataframes: list with dataframes with info on metabolites (list of dataframes)
get_list_dataframes_from_dir <- function(dir_with_subdirs) {
  # get a list of all metabolite files
  txt_files_paths <- list.files(dir_with_subdirs, pattern = "*.txt", recursive = FALSE, full.names = TRUE)
  # put all metabolites into one list
  list_of_dataframes <- lapply(txt_files_paths, read.table, sep = "\t", header = TRUE, quote = "")
  names(list_of_dataframes) <- gsub(".txt", "", basename(txt_files_paths))

  return(list_of_dataframes)
}

#' Merge patient Z-scores with metabolite info
#'
#' @param list_df_metabolite_groups: list of dataframes with metabolite information for different metabolite classes (list)
#' @param zscore_df: dataframe with metabolite Z-scores for all patient
#'
#' @return list_dfs_metabs_info_zscores: list of dataframes for each metabolite class
#' containing info and zscores for all samples
merge_metabolite_info_zscores <- function(list_df_metabolite_groups, zscore_df) {
  # remove HMDB_name column and "_Zscore" from column (patient) names
  zscore_df <- zscore_df %>%
    select(-HMDB_name)

  # put data into pages, max 20 violin plots per page in PDF
  list_dfs_metabs_info_zscores <- list()

  for (metabolite_class in names(list_df_metabolite_groups)) {
    # select the metabolite_class dataframe and select the HMDB_code and HMDB_name columns
    metabolite_info_df <- list_df_metabolite_groups[[metabolite_class]] %>% select(HMDB_code, HMDB_name)

    # Pad or truncate the HMDB names
    metabolite_info_df <- pad_truncate_hmdb_names(metabolite_info_df, 45, " ")

    # Join metabolite info with the Z-score dataframe
    metabolite_zscore_df <- metabolite_info_df %>%
      inner_join(zscore_df, by = "HMDB_code") %>%
      select(-HMDB_code)

    # put the data frame in long format
    metabolite_zscore_df_long <- reshape2::melt(
      metabolite_zscore_df,
      id.vars = "HMDB_name",
      variable.name = "Sample",
      value.name = "Z_score"
    )
    # Add the dataframe sorted on HMDB_name to a list
    list_dfs_metabs_info_zscores[[metabolite_class]] <- metabolite_zscore_df_long
  }

  return(list_dfs_metabs_info_zscores)
}

#' Combine patient and control data for each page of the violinplot pdf
#'
#' @param metab_interest_patients: list of dataframes with data for each metabolite and patient (list)
#' @param metab_interest_controls:  list of dataframes with data for each metabolite and control (list)
#' @param number_of_plots_per_page: number of plots per page in the violinplot pdf (integer)
#' @param number_of_patients: number of patients (integer)
#' @param number_of_controls: number of controls (integer)
#'
#' @return list_metabolite_df_per_page: list of dataframes with metabolite Z-scores for each patient and control,
#'          the length of list is the number of pages for the violinplot pdf (list)
get_data_per_metabolite_class <- function(
    metab_interest_patients,
    metab_interest_controls,
    number_of_plots_per_page,
    number_of_patients,
    number_of_controls) {
  list_metabolite_df_per_page <- list()
  metabolite_categories <- c()

  for (metabolite_class in names(metab_interest_patients)) {
    # Get the data for patients and controls for the metab_interest_sorted list
    metabolite_class_patients_df <- metab_interest_patients[[metabolite_class]]
    metabolite_class_controls_df <- metab_interest_controls[[metabolite_class]]

    # Get all metabolites and create list with HMDB names of max nr_plots_perpage long
    metabolites <- unique(metabolite_class_patients_df$HMDB_name)
    metabolites_in_chunks <- split(metabolites, ceiling(seq_along(metabolites) / number_of_plots_per_page))
    number_of_chunks_metabolites <- length(metabolites_in_chunks)

    # Get a list of plot data per page
    page_plot_data_list <- get_list_page_plot_data(
      metabolites_in_chunks,
      metabolite_class_patients_df,
      metabolite_class_controls_df,
      number_of_plots_per_page
    )

    # Add new items to main list
    list_metabolite_df_per_page <- append(list_metabolite_df_per_page, page_plot_data_list)
    # create list of page headers
    metabolite_categories <- c(metabolite_categories, paste(metabolite_class, seq(number_of_chunks_metabolites), sep = "_"))
  }
  # add page headers to list
  names(list_metabolite_df_per_page) <- metabolite_categories

  return(list_metabolite_df_per_page)
}

#' Get patient data to be uploaded to Helix
#'
#' @param list_dfs_metab_classes_zscores: list of dataframes with metabolite Z-scores for each sample/patient (list)
#' @param list_metabolite_classes: list of tables with metabolites for Helix and violin plots (list)
#'
#' @return df_zscores_to_helix: dataframe with patient data with only metabolites for Helix and violin plots
#'          with Helix name, high/low Z-score cutoffs
prepare_helix_patient_data <- function(list_dfs_metab_classes_zscores, list_metabolite_classes) {
  # Combine Z-scores of metab groups together
  metabolite_zscore_dataframe <- bind_rows(list_dfs_metab_classes_zscores)

  # Change the Sample column to characters, trim HMDB_name and split HMDB_name in new column
  metabolite_zscore_dataframe <- metabolite_zscore_dataframe %>%
    mutate(
      Sample = as.character(Sample),
      HMDB_name = str_trim(HMDB_name, "right"),
      HMDB_name_split = str_split_fixed(HMDB_name, "nitine;", 2)[, 1]
    )

  # Combine metabolite classes
  dims_helix_metabolite_df <- bind_rows(list_metabolite_classes)

  # Filter for Helix metabolites and split HMDB_name column for matching with metabolite_zscore_dataframe
  dims_helix_metabolite_df <- dims_helix_metabolite_df %>%
    filter(Helix == "ja") %>%
    mutate(HMDB_name_split = str_split_fixed(HMDB_name, "nitine;", 2)[, 1]) %>%
    select(HMDB_name_split, Helix_naam, high_zscore, low_zscore)

  # Filter DIMS results for metabolites for Helix and combine Helix info
  df_zscores_to_helix <- metabolite_zscore_dataframe %>%
    filter(HMDB_name_split %in% dims_helix_metabolite_df$HMDB_name_split) %>%
    left_join(dims_helix_metabolite_df, by = join_by(HMDB_name_split)) %>%
    select(HMDB_name, Sample, Z_score, Helix_naam, high_zscore, low_zscore)

  return(df_zscores_to_helix)
}

#' Getting the intensities for calculating ratio Z-scores
#' Retrieving a vector of intensities for a particular fraction side of the ratios for all samples.
#'
#' @param ratios_metabs_df: dataframe with HMDB codes for the ratios (dataframe)
#' @param row_index: index of the row in the ratios_metabs_df (integer)
#' @param intensities_zscore_df: dataframe with intensities for each sample (dataframe)
#' @param fraction_side: either numerator or denominator, which side of the fraction (string)
#' @param intensity_cols: names of the columns that contain the intensities (string)
#'
#' @returns fraction_side_intensity: a vector of intensities (vector of integers)
get_intensities_fraction_side <- function(ratios_metabs_df, row_index, intensities_zscore_df, fraction_side, intensity_cols) {
  # get the HMDB ID(s) for the given fraction side
  fraction_side_hmdb_ids <- ratios_metabs_df[row_index, fraction_side]
  if (grepl("plus", fraction_side_hmdb_ids)) {
    # if fraction side contains "plus", split to get both HMDB IDs
    fraction_side_hmdb_id_list <- strsplit(fraction_side_hmdb_ids, "plus")[[1]]
    # get intensities for both HMDB IDs for all samples
    fraction_side_intensity_list <- intensities_zscore_df %>%
      filter(HMDB_code %in% fraction_side_hmdb_id_list) %>%
      select(any_of(intensity_cols))
    # sum intensities to 1 intensity per samples
    fraction_side_intensity <- apply(fraction_side_intensity_list, 2, sum)
  } else if (fraction_side_hmdb_ids == "one") {
    # set intensity to 1
    fraction_side_intensity <- 1
  } else {
    # get intensities of the HMDB ID for all samples
    fraction_side_intensity <- intensities_zscore_df %>%
      filter(HMDB_code == fraction_side_hmdb_ids) %>%
      select(any_of(intensity_cols))
  }
  # vector of intensities for all samples
  fraction_side_intensity <- as.numeric(fraction_side_intensity)
  return(fraction_side_intensity)
}

#' Get the sample IDs for columns that have Z-score and intensities
#'
#' @param colnames_zscore_cols: vector of sample IDs from the dataframe containing Z-scores (vector of strings)
#' @param colnames_intensity_cols: vector of sample IDs form the dataframe containing intensities (vector of strings)
#'
#' @returns colnames_intersect: vector of sample IDs that are in both input vectors, ending on "_Zscore" (vector of strings)
get_sample_ids_with_zscores <- function(colnames_zscore_cols, colnames_intensity_cols) {
  colnames_intersect <- intersect(
    paste0(colnames_intensity_cols, "_Zscore"),
    grep("_Zscore", colnames_zscore_cols, value = TRUE)
  )
  return(colnames_intersect)
}

#' Pad or truncate HMDB names to a fixed width
#' Add spaces or remove HMDB name characters till the length of the name equals the 'width'
#'
#' @param metabolite_info_df: A dataframe containing a column `HMDB_name` (character).
#' @param width: Integer target width for the display names. Default is 45.
#' @param pad_character: Single character used for padding. Default is a space `" "`.
#'
#' @return metabolite_info_df: A dataframe where the HMDB names are transformed
pad_truncate_hmdb_names <- function(metabolite_info_df, width, pad_character) {
  # Change the HMDB_name column so all names have 45 characters
  # remove characters if name is longer and add "..."
  # add empty spaces till 45 charachters if name is shorter
  # keep the name if name is exactly 45 characters
  keep_lenght <- width - 3
  metabolite_info_df <- metabolite_info_df %>% mutate(HMDB_name = case_when(
    str_length(HMDB_name) > width ~ str_c(str_sub(HMDB_name, 1, keep_lenght), "..."),
    str_length(HMDB_name) < width ~ str_pad(HMDB_name, width, side = "right", pad = pad_character),
    TRUE ~ HMDB_name
  ))
  return(metabolite_info_df)
}

#' Get a list of dataframes for each chunk
#' For each chunk, get a dataframe containing the metabolites in that chunk and add it to the list
#'
#' @param metabolites_in_chunks: list of vectors, each containing metabolites
#' @param metabolite_class_patients_df: dataframe of Z-scores for all patient
#' @param metabolite_class_controls_df: dataframe of Z-scores for all control
#' @param number_of_plots_per_page: integer containing the number of metabolites per plot per page
#'
#' @returns page_plot_data_list: a list of dataframes containing Z-scores
get_list_page_plot_data <- function(
    metabolites_in_chunks,
    metabolite_class_patients_df,
    metabolite_class_controls_df,
    number_of_plots_per_page) {
  # For each chunk, get a dataframe containing the metabolites and add to a list
  page_plot_data_list <- lapply(metabolites_in_chunks, function(metabolite_names_chunk) {
    patients_df_chunk <- metabolite_class_patients_df %>% filter(HMDB_name %in% metabolite_names_chunk)
    controls_df_chunk <- metabolite_class_controls_df %>% filter(HMDB_name %in% metabolite_names_chunk)

    # Combine both dataframes
    patients_controls_df_chunk <- rbind(patients_df_chunk, controls_df_chunk)
    metabolite_order <- make_metabolite_order(number_of_plots_per_page, metabolite_names_chunk)

    # Set the order of the metabolites for the violin plots
    attr(patients_controls_df_chunk, "y_order") <- rev(metabolite_order)

    return(patients_controls_df_chunk)
  })
}

#' Make the order of metabolites for the violin plots
#' Create the order of metabolites and add empty strings if the number of metabolites is lower than
#' the number of plots per page.
#'
#' @param number_of_plots_per_page: integer containing the number of metabolites per plot per page
#' @param metabolite_names_chunk: list of vectors, each containing metabolites
#'
#' @returns metabolite_order: a vector containing all metabolites and possibly empty strings
make_metabolite_order <- function(number_of_plots_per_page, metabolite_names_chunk) {
  # Add empty dummy's to extend the number of metabs to the nr_plots_perpage
  number_of_plots_missing <- number_of_plots_per_page - length(metabolite_names_chunk)
  if (number_of_plots_missing > 0) {
    dummy_names <- paste0(" ", strrep(" ", seq_len(number_of_plots_missing)))
    metabolite_order <- c(metabolite_names_chunk, dummy_names)
  } else {
    metabolite_order <- metabolite_names_chunk
  }
  return(metabolite_order)
}

#' Check for Diagnostics patients with correct patient number (e.g. starting with "P2024M")
#'
#' @param patient_column: a column from dataframe with IDs (character vector)
#'
#' @return: a logical vector with TRUE or FALSE for each element (vector)
is_diagnostic_patients <- function(patient_column) {
  diagnostic_patients <- grepl("^P[0-9]{4}M", patient_column)

  return(diagnostic_patients)
}

#' Get the output dataframe for Helix
#'
#' @param protocol_name: protocol name (string)
#' @param df_metabs_helix: dataframe with metabolite Z-scores for patients (dataframe)
#'
#' @return: dataframe with patient metabolite Z-scores in correct format for Helix
transform_metab_df_to_helix_df <- function(protocol_name, df_metabs_helix) {
  # Remove positive controls
  df_metabs_helix <- df_metabs_helix %>% filter(is_diagnostic_patients(Sample))

  # Add 'Vial' column, each patient has unique ID
  df_metabs_helix <- df_metabs_helix %>%
    group_by(Sample) %>%
    mutate(Vial = cur_group_id()) %>%
    ungroup()

  # Split patient number into labnummer and Onderzoeksnummer
  df_metabs_helix <- add_lab_id_and_onderzoeksnr(df_metabs_helix)

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
add_lab_id_and_onderzoeksnr <- function(df_metabs_helix) {
  # Split patient number into labnummer and Onderzoeksnummer
  for (row in seq_len(nrow(df_metabs_helix))) {
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
get_top_metabolites_df <- function(patient_name, dims_helix_table) {
  # extract data for patient of interest (patient_name)
  patient_metabs_helix <- dims_helix_table %>%
    filter(Sample == patient_name) %>%
    mutate(Z_score = round(Z_score, 2))

  patient_high_df <- patient_metabs_helix %>% filter(Z_score > high_zscore)
  patient_low_df <- patient_metabs_helix %>% filter(Z_score < low_zscore)

  if (nrow(patient_high_df) > 0 || nrow(patient_low_df) > 0) {
    # sort tables on zscore
    patient_high_df <- patient_high_df %>%
      arrange(desc(Z_score)) %>%
      select(c(HMDB_name, Z_score))
    patient_low_df <- patient_low_df %>%
      arrange(Z_score) %>%
      select(c(HMDB_name, Z_score))
  }
  # add lines for increased, decreased
  line_increased <- c("Increased", "")
  line_decreased <- c("Decreased", "")

  # combine the two lists
  top_metab_patient <- rbind(line_increased, patient_high_df, line_decreased, patient_low_df)

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
prepare_toplist <- function(patient_id, zscore_patients, num_of_highest_metabolites, num_of_lowest_metabolites) {
  patient_df <- zscore_patients %>%
    select(HMDB_code, HMDB_name, !!sym(patient_id)) %>%
    arrange(!!sym(patient_id))

  # Get lowest Zscores
  patient_df_low <- patient_df[1:num_of_lowest_metabolites, ]
  patient_df_low <- patient_df_low %>% mutate(across(!!sym(patient_id), ~ round(.x, 2)))

  # Get highest Zscores
  patient_df_high <- patient_df[nrow(patient_df):(nrow(patient_df) - num_of_highest_metabolites + 1), ]
  patient_df_high <- patient_df_high %>% mutate(across(!!sym(patient_id), ~ round(.x, 2)))

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

  # get the names and numbers in the table aligned
  table_theme <- ttheme_default(
    core = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 6)),
    colhead = list(fg_params = list(fontsize = 8, fontface = "bold"))
  )

  # patient plots, create the PDF device
  patient_id_sub <- patient_id
  suffix <- ""
  if (grepl("Diagnostics", pdf_dir) && is_diagnostic_patients(patient_id)) {
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
    height = plot_height
  )

  # page headers:
  page_headers <- names(metab_perpage)

  # put table into PDF file, if not empty
  if (!is.null(dim(top_metab_pt))) {
    max_rows_per_page <- 35
    total_rows <- nrow(top_metab_pt)
    number_of_pages <- ceiling(total_rows / max_rows_per_page)

    # get the names and numbers in the table aligned
    table_theme <- ttheme_default(
      core = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 6)),
      colhead = list(fg_params = list(fontsize = 8, fontface = "bold"))
    )

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
    geom_point(
      data = patient_zscore_df, aes(color = Z_score),
      size = 3.5 * circlesize, shape = 22, fill = "white", na.rm = TRUE
    ) +
    # Add the Z-score at the right side of the plot
    geom_text(
      data = patient_zscore_df,
      aes(16, label = paste0("Z=", round(Z_score, 2))),
      hjust = "left", vjust = +0.2, size = 3, na.rm = TRUE
    ) +
    # Set colour for the Z-score of the selected patient
    scale_fill_gradientn(
      colors = colors_plot, values = NULL, space = "Lab",
      na.value = "grey50", guide = "colourbar", aesthetics = "colour"
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
#' @param expected_biomarkers_df: dataframe with information for HMDB codes about IEMs (dataframe)
#' @param zscore_patients: dataframe containing Z-scores for patient (dataframe)
#' @param sample_cols: vector containing column names with intensities and Z-scores for patients (vector)
#'
#' @returns probability_score: a dataframe with probability scores for IEMs for each patient (dataframe)
run_diem_algorithm <- function(expected_biomarkers_df, zscore_patients_df, sample_cols) {
  # Rank the metabolites for each patient individually
  ranking_patients <- zscore_patients_df %>%
    mutate(across(-c(HMDB_code, HMDB_name), rank_patient_zscores))

  ranking_patients <- merge(
    x = expected_biomarkers_df, y = ranking_patients,
    by.x = c("HMDB_code"), by.y = c("HMDB_code")
  )

  zscore_expected_df <- merge(
    x = expected_biomarkers_df, y = zscore_patients_df,
    by.x = c("HMDB_code"), by.y = c("HMDB_code")
  )

  # Reduce sample_cols to contain only the patient names
  if (sum(grepl("Zscore", sample_cols)) > 0) {
    sample_cols <- sample_cols[-grep("Zscore", sample_cols)]
  }

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
    filter(!duplicated(select(., Disease, M.z)) |
             !duplicated(select(., Disease, M.z), fromLast = FALSE)) %>%
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
#' @param run_name: name of the run, for the file name (string)
save_prob_scores_to_excel <- function(probability_score, run_name) {
  # Create conditional formatting for output Excel sheet. Colors according to values.
  wb <- createWorkbook()
  addWorksheet(wb, "Probability Scores")
  writeData(wb, "Probability Scores", probability_score)
  conditionalFormatting(wb, "Probability Scores",
    cols = 2:ncol(probability_score), rows = seq_len(nrow(probability_score)),
    type = "colourScale", style = c("white", "#FFFDA2", "red"), rule = c(1, 10, 100)
  )
  saveWorkbook(wb, file = paste0("./dIEM_algoritme_output_", run_name, ".xlsx"), overwrite = TRUE)
  rm(wb)
}

#' Make and save dIEM plots
#'
#' @param diem_probability_score: dataframe with dIEM probability scores
#' @param patient_col_names: vector containing all patient column names
#' @param expected_biomarkers_df: dataframe with information for HMDB codes about IEMs
#' @param zscore_patients_df: dataframe containing Z-scores for all patients
#' @param zscore_controls_df: dataframe containing Z-scores for all controls
#' @param nr_plots_perpage: integer containing the number of metabolites per page
#' @param number_of_samples: list containing the number of patients and controls
#' @param number_of_metabolites: list containing the number of metabolites for the top and lowest table
#'
#' @returns patient_no_iem: vector of patient IDs that have no IEMs
make_and_save_diem_plots <- function(
    diem_probability_score,
    patient_col_names,
    expected_biomarkers_df,
    zscore_patients_df,
    zscore_controls_df,
    nr_plots_perpage,
    number_of_samples,
    number_of_metabolites,
    iem_variables,
    explanation_violin_plot) {
  diem_plot_dir <- paste("./dIEM_plots", sep = "/")
  dir.create(diem_plot_dir)

  # reduce patient_col_names to contain only the patient names
  if (sum(grepl("Zscore", patient_col_names)) > 0) {
    patient_col_names <- patient_col_names[-grep("Zscore", patient_col_names)]
  }

  patient_no_iem <- c()

  for (patient_id in patient_col_names) {
    # Select the top IEMs and filter on the IEM threshold
    patient_top_iems_probs <- diem_probability_score %>%
      select(c(Disease, !!sym(patient_id))) %>%
      arrange(desc(!!sym(patient_id))) %>%
      slice(1:iem_variables$top_number_iem_diseases) %>%
      filter(!!sym(patient_id) >= iem_variables$threshold_iem)

    if (nrow(patient_top_iems_probs) > 0) {
      list_metabolites_top_iems <- get_probabilities_top_iems(patient_top_iems_probs, expected_biomarkers_df, patient_id)

      # Get the Z-scores with metabolite information
      metabolites_iem_sorted <- merge_metabolite_info_zscores(list_metabolites_top_iems, zscore_patients_df)
      metabolites_iem_controls <- merge_metabolite_info_zscores(list_metabolites_top_iems, zscore_controls_df)

      # Get a list of dataframes for each IEM
      diem_metabolites_perpage <- get_data_per_metabolite_class(
        metabolites_iem_sorted,
        metabolites_iem_controls,
        nr_plots_perpage,
        number_of_samples$patients,
        number_of_samples$controls
      )
      # Get a dataframe of the top metabolites
      top_metabolites_patient <- prepare_toplist(
        patient_id,
        zscore_patients_df,
        number_of_metabolites$highest,
        number_of_metabolites$lowest
      )

      # Generate and save dIEM violin plots
      create_pdf_violin_plots(
        diem_plot_dir,
        patient_id,
        diem_metabolites_perpage,
        top_metabolites_patient,
        explanation_violin_plot
      )
    } else {
      patient_no_iem <- c(patient_no_iem, patient_id)
    }
  }
  return(patient_no_iem)
}

#' Get the IEM probabilities for a patient for all diseases
#'
#' @param patient_top_iems_probs: dataframe containing the probability scores for diseases for a patient
#' @param expected_biomarkers_df: dataframe with information for HMDB codes about IEMs
#' @param patient_id: string containing the patien ID
#'
#' @returns list_metabolites_iems: list of dataframes containing the HMDB codes and names for all diseases
get_probabilities_top_iems <- function(patient_top_iems_probs, expected_biomarkers_df, patient_id) {
  # Get the metabolites for each IEM and their probability
  list_metabolites_iems <- list()
  metabolites_iems_names <- c()

  for (iem in patient_top_iems_probs$Disease) {
    # get the IEM probabilities for the selected patient
    iem_probablity <- patient_top_iems_probs %>%
      filter(Disease == iem) %>%
      pull(!!sym(patient_id))
    metabolites_iems_names <- c(metabolites_iems_names, paste0(iem, ", probability score ", iem_probablity))
    # get the HMDB codes and names for the IEM of the selected patient
    metabolites_iem_df <- expected_biomarkers_df %>%
      filter(Disease == iem) %>%
      select(HMDB_code, HMDB_name)
    # Add for each IEM the HMDB codes and names to a list
    list_metabolites_iems[[iem]] <- metabolites_iem_df
  }
  names(list_metabolites_iems) <- metabolites_iems_names

  return(list_metabolites_iems)
}

#' Save a list of patient IDs to a text file
#'
#' @param threshold_iem: integer containing the IEM threshold
#' @param patient_no_iem: vector containing patient IDs
save_patient_no_iem <- function(threshold_iem, patient_no_iem) {
  patient_no_iem <- c(
    paste0(
      "The following patient(s) did not have dIEM probability scores higher than ",
      threshold_iem, " :"
    ),
    patient_no_iem
  )
  write(file = paste0("./missing_probability_scores.txt"), patient_no_iem)
}
