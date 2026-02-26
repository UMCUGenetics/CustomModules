# unit tests for GenerateViolinPlots

suppressPackageStartupMessages(library("dplyr"))
library(reshape2)
library(openxlsx)
library(ggplot2)
suppressPackageStartupMessages(library("gridExtra"))
library(stringr)

source("../../export/generate_violin_plots_functions.R")

testthat::test_that("get_intensities_fraction_side: Get intensities for calculating the ratios", {
  test_intensities_zscore_df <- read.delim(test_path("fixtures", "test_intensities_zscore_df.txt"))
  test_intensity_cols <- c("C101.1", "C102.1", "C103.1", "C104.1", "C105.1",
                           "P2025M1", "P2025M2", "P2025M3", "P2025M4", "P2025M5")
  test_ratios_metabs_df <- data.frame(
    HMDB.code = c("HMDBAA1", "HMDBAA2", "HMDBAB1"),
    Ratio_name = c("ratio1", "ratio2", "ratio3"),
    HMDB_numerator = c("HMDB001", "HMDB002plusHMDB003", "HMDB004"),
    HMDB_denominator = c("HMDB011", "HMDB012", "one")
  )
  expect_equal(
    get_intensities_fraction_side(
      test_ratios_metabs_df,
      2,
      test_intensities_zscore_df,
      "HMDB_numerator",
      test_intensity_cols
    ),
    c(2500, 2625, 1750, 2525, 2700, 2500, 2625, 1750, 2525, 2700)
  )
  expect_equal(
    get_intensities_fraction_side(
      test_ratios_metabs_df,
      3,
      test_intensities_zscore_df,
      "HMDB_denominator",
      test_intensity_cols
    ),
    1
  )
})

testthat::test_that("get_sample_ids_with_zscores: Get samples that have both intensity and Z-score columns", {
  test_intensities_zscore_df <- read.delim(test_path("fixtures", "test_intensities_zscore_df.txt"))

  test_intensity_cols <- c(
    "C101.1", "C102.1", "C103.1", "C104.1", "C105.1",
    "P2025M1", "P2025M2", "P2025M3", "P2025M4", "P2025M5"
  )

  expect_equal(
    length(get_sample_ids_with_zscores(colnames(test_intensities_zscore_df), test_intensity_cols)),
    10
  )
  expect_equal(
    get_sample_ids_with_zscores(colnames(test_intensities_zscore_df), test_intensity_cols),
    c(
      "C101.1_Zscore", "C102.1_Zscore", "C103.1_Zscore", "C104.1_Zscore", "C105.1_Zscore",
      "P2025M1_Zscore", "P2025M2_Zscore", "P2025M3_Zscore", "P2025M4_Zscore", "P2025M5_Zscore"
    )
  )

  test_intensity_cols <- c("C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "P2025M1", "P2025M2", "P2025M3")
  expect_equal(
    length(get_sample_ids_with_zscores(colnames(test_intensities_zscore_df), test_intensity_cols)),
    8
  )
  expect_equal(
    get_sample_ids_with_zscores(colnames(test_intensities_zscore_df), test_intensity_cols),
    c(
      "C101.1_Zscore", "C102.1_Zscore", "C103.1_Zscore", "C104.1_Zscore", "C105.1_Zscore",
      "P2025M1_Zscore", "P2025M2_Zscore", "P2025M3_Zscore"
    )
  )
})

testthat::test_that("get_list_dataframes_from_dir: Get a list with dataframes from metabilite files in a directory", {
  test_path_metabolite_groups <- test_path("fixtures/test_metabolite_groups/Diagnostics")

  expect_type(get_list_dataframes_from_dir(test_path_metabolite_groups), "list")
  expect_identical(
    names(get_list_dataframes_from_dir(test_path_metabolite_groups)),
    c("test_acyl_carnitines", "test_crea_gua")
  )

  expect_identical(
    colnames(get_list_dataframes_from_dir(test_path_metabolite_groups)$test_acyl_carnitines),
    c("HMDB_code", "HMDB_name", "Helix", "Helix_naam", "high_zscore", "low_zscore")
  )
  expect_identical(
    get_list_dataframes_from_dir(test_path_metabolite_groups)$test_acyl_carnitines$HMDB_name,
    c("metab1", "metab3", "ratio1")
  )
  expect_identical(
    get_list_dataframes_from_dir(test_path_metabolite_groups)$test_acyl_carnitines$Helix,
    c("ja", "nee", "ja")
  )

  expect_identical(
    colnames(get_list_dataframes_from_dir(test_path_metabolite_groups)$test_crea_gua),
    c("HMDB_code", "HMDB_name", "Helix", "Helix_naam", "high_zscore", "low_zscore")
  )
  expect_identical(
    get_list_dataframes_from_dir(test_path_metabolite_groups)$test_crea_gua$HMDB_name,
    c("metab4", "metab11")
  )
  expect_identical(
    get_list_dataframes_from_dir(test_path_metabolite_groups)$test_crea_gua$Helix,
    c("ja", "ja")
  )
})

testthat::test_that("merge_metabolite_info_zscores: Combine metabolite info dataframe and Z-score dataframe", {
  test_acyl_carnitines_df <- read.delim(test_path("fixtures/test_metabolite_groups/Diagnostics", "test_acyl_carnitines.txt"))
  test_crea_gua_df <- read.delim(test_path("fixtures/test_metabolite_groups/Diagnostics", "test_crea_gua.txt"))

  test_metab_list_all <- list(test_acyl_carnitines_df, test_crea_gua_df)
  names(test_metab_list_all) <- c("test_acyl_carnitines", "test_crea_gua")

  test_patient_cols <- c("P2025M1_Zscore", "P2025M2_Zscore", "P2025M3_Zscore", "P2025M4_Zscore", "P2025M5_Zscore")
  test_intensities_zscore_df <- read.delim(test_path("fixtures", "test_intensities_zscore_df.txt"))
  test_zscore_patients_df <- test_intensities_zscore_df %>% select(HMDB_code, HMDB_name, any_of(test_patient_cols))

  expect_type(merge_metabolite_info_zscores(test_metab_list_all, test_zscore_patients_df), "list")
  expect_identical(
    names(merge_metabolite_info_zscores(test_metab_list_all, test_zscore_patients_df)),
    c("test_acyl_carnitines", "test_crea_gua")
  )

  expect_identical(
    colnames(merge_metabolite_info_zscores(test_metab_list_all, test_zscore_patients_df)$test_acyl_carnitines),
    c("HMDB_name", "Sample", "Z_score")
  )
  expect_identical(
    (merge_metabolite_info_zscores(
      test_metab_list_all,
      test_zscore_patients_df
    )$test_acyl_carnitines$Z_score),
    c(0.31, 2.34, 2.45, 1.45, 2.14, -1.44, 12.18, -0.18, 3.22, -3.18)
  )
  expect_identical(
    as.character(merge_metabolite_info_zscores(
      test_metab_list_all,
      test_zscore_patients_df
    )$test_acyl_carnitines$Sample),
    c(
      "P2025M1_Zscore", "P2025M1_Zscore", "P2025M2_Zscore", "P2025M2_Zscore", "P2025M3_Zscore",
      "P2025M3_Zscore", "P2025M4_Zscore", "P2025M4_Zscore", "P2025M5_Zscore", "P2025M5_Zscore"
    )
  )
  expect_equal(
    nchar(merge_metabolite_info_zscores(
      test_metab_list_all,
      test_zscore_patients_df
    )$test_acyl_carnitines$HMDB_name[1]),
    45
  )

  expect_identical(
    colnames(merge_metabolite_info_zscores(test_metab_list_all, test_zscore_patients_df)$test_crea_gua),
    c("HMDB_name", "Sample", "Z_score")
  )
  expect_identical(
    (merge_metabolite_info_zscores(
      test_metab_list_all,
      test_zscore_patients_df
    )$test_crea_gua$Z_score),
    c(0.84, -0.46, -0.15, -1.51, -0.78, 1.68, 0.84, 1.48, 0.47, 1.18)
  )
  expect_identical(
    as.character(merge_metabolite_info_zscores(
      test_metab_list_all,
      test_zscore_patients_df
    )$test_crea_gua$Sample),
    c(
      "P2025M1_Zscore", "P2025M1_Zscore", "P2025M2_Zscore", "P2025M2_Zscore", "P2025M3_Zscore",
      "P2025M3_Zscore", "P2025M4_Zscore", "P2025M4_Zscore", "P2025M5_Zscore", "P2025M5_Zscore"
    )
  )
  expect_equal(
    nchar(merge_metabolite_info_zscores(
      test_metab_list_all,
      test_zscore_patients_df
    )$test_crea_gua$HMDB_name[1]),
    45
  )
})

testthat::test_that("get_data_per_metabolite_class: Combine patient and control data for each page of the violinplot pdf", {
  test_acyl_carnitines_pat <- read.delim(test_path("fixtures/", "test_acyl_carnitines_patients.txt"))
  test_acyl_carnitines_ctrl <- read.delim(test_path("fixtures/", "test_acyl_carnitines_controls.txt"))

  test_crea_gua_pat <- read.delim(test_path("fixtures/", "test_crea_gua_patients.txt"))
  test_crea_gua_ctrl <- read.delim(test_path("fixtures/", "test_crea_gua_controls.txt"))

  test_metab_interest_sorted <- list(test_acyl_carnitines_pat, test_crea_gua_pat)
  names(test_metab_interest_sorted) <- c("test_acyl_carnitines", "test_crea_gua")

  test_metab_interest_contr <- list(test_acyl_carnitines_ctrl, test_crea_gua_ctrl)
  names(test_metab_interest_contr) <- c("test_acyl_carnitines", "test_crea_gua")

  test_nr_plots_perpage <- 1
  test_nr_pat <- 5
  test_nr_contr <- 5

  expect_type(
    get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat, test_nr_contr
    ),
    "list"
  )
  expect_equal(
    length(get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat, test_nr_contr
    )),
    4
  )
  expect_identical(
    names(get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat, test_nr_contr
    )),
    c("test_acyl_carnitines_1", "test_acyl_carnitines_2", "test_crea_gua_1", "test_crea_gua_2")
  )
  expect_identical(
    unique(get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat,
      test_nr_contr
    )$test_acyl_carnitines_1$HMDB_name),
    c("metab1                                       ")
  )
  expect_identical(
    get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat,
      test_nr_contr
    )$test_acyl_carnitines_1$Sample,
    c("P2025M1", "P2025M2", "P2025M3", "P2025M4", "P2025M5", "C101.1", "C102.1", "C103.1", "C104.1", "C105.1")
  )

  test_nr_plots_perpage <- 2

  expect_equal(
    length(get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat, test_nr_contr
    )),
    2
  )
  expect_identical(
    names(get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat, test_nr_contr
    )),
    c("test_acyl_carnitines_1", "test_crea_gua_1")
  )
  expect_identical(
    unique(get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat,
      test_nr_contr
    )$test_acyl_carnitines_1$HMDB_name),
    c("metab1                                       ", "metab3                                       ")
  )
  expect_identical(
    get_data_per_metabolite_class(
      test_metab_interest_sorted, test_metab_interest_contr,
      test_nr_plots_perpage, test_nr_pat,
      test_nr_contr
    )$test_acyl_carnitines_1$Sample,
    c(
      "P2025M1", "P2025M1", "P2025M2", "P2025M2", "P2025M3",
      "P2025M3", "P2025M4", "P2025M4", "P2025M5", "P2025M5",
      "C101.1", "C101.1", "C102.1", "C102.1", "C103.1", "C103.1", "C104.1", "C104.1", "C105.1", "C105.1"
    )
  )
})

testthat::test_that("prepare_helix_patient_data: Generate a dataframe with information for Helix", {
  test_acyl_carnitines_pat <- read.delim(test_path("fixtures", "test_acyl_carnitines_patients.txt"))
  test_crea_gua_pat <- read.delim(test_path("fixtures", "test_crea_gua_patients.txt"))

  test_metab_interest_sorted <- list(test_acyl_carnitines_pat, test_crea_gua_pat)
  names(test_metab_interest_sorted) <- c("test_acyl_carnitines", "test_crea_gua")

  test_acyl_carnitines_df <- read.delim(test_path("fixtures/test_metabolite_groups/Diagnostics", "test_acyl_carnitines.txt"))
  test_crea_gua_df <- read.delim(test_path("fixtures/test_metabolite_groups/Diagnostics", "test_crea_gua.txt"))

  test_metab_list_all <- list(test_acyl_carnitines_df, test_crea_gua_df)
  names(test_metab_list_all) <- c("test_acyl_carnitines", "test_crea_gua")

  expect_identical(
    colnames(prepare_helix_patient_data(test_metab_interest_sorted, test_metab_list_all)),
    c("HMDB_name", "Sample", "Z_score", "Helix_naam", "high_zscore", "low_zscore")
  )
  expect_equal(
    dim(prepare_helix_patient_data(test_metab_interest_sorted, test_metab_list_all)),
    c(15, 6)
  )
  expect_identical(
    unique(prepare_helix_patient_data(test_metab_interest_sorted, test_metab_list_all)$HMDB_name),
    c("metab1", "metab4", "metab11")
  )
  expect_false("ratio1" %in% prepare_helix_patient_data(test_metab_interest_sorted, test_metab_list_all)$HMDB_name)
  expect_equal(
    prepare_helix_patient_data(test_metab_interest_sorted, test_metab_list_all)$Z_score,
    c(0.31, 2.45, 2.14, 12.18, 3.22, 0.84, -0.46, -0.15, -1.51, -0.78, 1.68, 0.84, 1.48, 0.47, 1.18)
  )
})

testthat::test_that("is_diagnostic_patients: Check for diagnostic patients", {
  test_patient_column <- c("P2025M1", "P2025M2", "P2025M3", "P2025M4", "C101.1", "C102.1", "P2025D1", "P225M1")

  expect_equal(is_diagnostic_patients(test_patient_column), c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(length(is_diagnostic_patients(test_patient_column)), 8)
})

testthat::test_that("add_lab_id_and_onderzoeksnr:Adding labnummer and Onderzoeksnummer to the Helix dataframe", {
  test_df_metabs_helix <- read.delim(test_path("fixtures/", "test_df_metabs_helix.txt"))

  test_df_metabs_helix <- test_df_metabs_helix %>%
    group_by(Sample) %>%
    mutate(Vial = cur_group_id()) %>%
    ungroup()

  expect_true("labnummer" %in% colnames(add_lab_id_and_onderzoeksnr(test_df_metabs_helix)))
  expect_true("Onderzoeksnummer" %in% colnames(add_lab_id_and_onderzoeksnr(test_df_metabs_helix)))
  expect_identical(
    unique(add_lab_id_and_onderzoeksnr(test_df_metabs_helix)$labnummer),
    c("2025M1", "2025M2", "2025M3", "2025M4", "2025M5")
  )
  expect_identical(
    unique(add_lab_id_and_onderzoeksnr(test_df_metabs_helix)$Onderzoeksnummer),
    c("MB2025/1", "MB2025/2", "MB2025/3", "MB2025/4", "MB2025/5")
  )
})

testthat::test_that("transform_metab_df_to_helix_df: Make the output for Helix", {
  test_protocol_name <- "test_protocol_name"

  test_df_metabs_helix <- read.delim(test_path("fixtures/", "test_df_metabs_helix.txt"))

  expect_equal(
    dim(transform_metab_df_to_helix_df(test_protocol_name, test_df_metabs_helix)),
    c(15, 6)
  )
  expect_identical(
    colnames(transform_metab_df_to_helix_df(test_protocol_name, test_df_metabs_helix)),
    c("Vial", "labnummer", "Onderzoeksnummer", "Protocol", "Name", "Amount")
  )
  expect_identical(
    unique(transform_metab_df_to_helix_df(test_protocol_name, test_df_metabs_helix)$Protocol),
    "test_protocol_name"
  )
  expect_identical(
    unique(transform_metab_df_to_helix_df(test_protocol_name, test_df_metabs_helix)$labnummer),
    c("2025M1", "2025M2", "2025M3", "2025M4", "2025M5")
  )
  expect_identical(
    unique(transform_metab_df_to_helix_df(test_protocol_name, test_df_metabs_helix)$Onderzoeksnummer),
    c("MB2025/1", "MB2025/2", "MB2025/3", "MB2025/4", "MB2025/5")
  )
  expect_equal(
    transform_metab_df_to_helix_df(test_protocol_name, test_df_metabs_helix)$Amount,
    c(0.31, 2.45, 2.14, 12.18, 3.22, 0.84, -0.46, -0.15, -1.51, -0.78, 1.68, 0.84, 1.48, 0.47, 1.18)
  )
})

testthat::test_that("get_top_metabolites_df: Create a dataframe with all metabolites that exceed the min and max Z-score cutoff", {
  test_df_metabs_helix <- read.delim(test_path("fixtures/", "test_df_metabs_helix.txt"))
  test_patient_id <- "P2025M1"

  expect_equal(
    dim(get_top_metabolites_df(test_patient_id, test_df_metabs_helix)),
    c(2, 2)
  )
  expect_equal(
    colnames(get_top_metabolites_df(test_patient_id, test_df_metabs_helix)),
    c("Metabolite", "Z-score")
  )

  test_patient_id <- "P2025M2"
  expect_equal(
    dim(get_top_metabolites_df(test_patient_id, test_df_metabs_helix)),
    c(4, 2)
  )
  expect_equal(
    get_top_metabolites_df(test_patient_id, test_df_metabs_helix)$`Z-score`,
    c("", "2.45", "", "-1.51")
  )

  test_patient_id <- "P2025M4"
  expect_equal(
    dim(get_top_metabolites_df(test_patient_id, test_df_metabs_helix)),
    c(3, 2)
  )
  expect_equal(
    get_top_metabolites_df(test_patient_id, test_df_metabs_helix)$`Z-score`,
    c("", "12.18", "")
  )
})

testthat::test_that("prepare_toplist: Create a dataframe with the top 20 highest and top 10 lowest metabolites", {
  test_zscore_patient_df <- read.delim(test_path("fixtures/", "test_zscore_patient_df.txt"))
  test_patient_id <- "P2025M1"
  test_num_of_highest_metabolites <- 20
  test_num_of_lowest_metabolites <- 10

  expect_equal(
    dim(prepare_toplist(
      test_patient_id,
      test_zscore_patient_df,
      test_num_of_highest_metabolites,
      test_num_of_lowest_metabolites
    )),
    c(32, 3)
  )
  expect_equal(
    colnames(prepare_toplist(
      test_patient_id,
      test_zscore_patient_df,
      test_num_of_highest_metabolites,
      test_num_of_lowest_metabolites
    )),
    c("HMDB_ID", "Metabolite", "Z-score")
  )
  expect_equal(
    prepare_toplist(
      test_patient_id,
      test_zscore_patient_df,
      test_num_of_highest_metabolites,
      test_num_of_lowest_metabolites
    )$HMDB_ID,
    c(
      "Increased", "HMDB030", "HMDB029", "HMDB028", "HMDB027", "HMDB026", "HMDB025", "HMDB024", "HMDB023",
      "HMDB022", "HMDB021", "HMDB020", "HMDB019", "HMDB018", "HMDB017", "HMDB016", "HMDB015", "HMDB014", "HMDB013",
      "HMDB012", "HMDB011", "Decreased", "HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB005", "HMDB006",
      "HMDB007", "HMDB008", "HMDB009", "HMDB010"
    )
  )
  expect_equal(
    prepare_toplist(
      test_patient_id,
      test_zscore_patient_df,
      test_num_of_highest_metabolites,
      test_num_of_lowest_metabolites
    )$`Z-score`,
    c(
      "", "30", "29", "28", "27", "26", "25", "24", "23", "22", "21", "20", "19", "18", "17", "16", "15",
      "14", "13", "12", "11", "", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"
    )
  )

  test_patient_id <- "P2025M2"

  expect_equal(
    prepare_toplist(
      test_patient_id,
      test_zscore_patient_df,
      test_num_of_highest_metabolites,
      test_num_of_lowest_metabolites
    )$Metabolite,
    c(
      "", "metab1", "metab2", "metab3", "metab4", "metab5", "metab6", "metab7", "metab8", "metab9", "metab10",
      "metab11", "metab12", "metab13", "metab14", "metab15", "metab16", "metab17", "metab18", "metab19", "metab20",
      "", "metab30", "metab29", "metab28", "metab27", "metab26", "metab25", "metab24", "metab23", "metab22",
      "metab21"
    )
  )
  expect_equal(
    prepare_toplist(
      test_patient_id,
      test_zscore_patient_df,
      test_num_of_highest_metabolites,
      test_num_of_lowest_metabolites
    )$`Z-score`,
    c(
      "", "-1", "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9", "-10", "-11", "-12", "-13", "-14", "-15",
      "-16", "-17", "-18", "-19", "-20", "", "-30", "-29", "-28", "-27", "-26", "-25", "-24", "-23", "-22", "-21"
    )
  )
})

testthat::test_that("create_pdf_violin_plots: Create a pdf with a table of top metabolites and violin plots", {
  local_edition(3)
  temp_dir <- "./"
  dir.create(paste0(temp_dir, "violin_plots/"))

  test_pdf_dir <- paste0(temp_dir, "violin_plots/")
  test_patient_id <- "P2025M1"
  test_explanation <- "Unit test Generate Violin Plots"

  test_acyl_carnitines_df <- read.delim(test_path("fixtures/", "test_acyl_carnitines_df.txt"))
  attr(test_acyl_carnitines_df, "y_order") <- rev(unique(test_acyl_carnitines_df$HMDB_name))
  test_crea_gua_df <- read.delim(test_path("fixtures/", "test_crea_gua_df.txt"))
  attr(test_crea_gua_df, "y_order") <- rev(unique(test_crea_gua_df$HMDB_name))

  test_metab_perpage <- list(test_acyl_carnitines_df, test_crea_gua_df)
  names(test_metab_perpage) <- c("test_acyl_carnitines", "test_crea_gua")

  test_top_metab_pt <- data.frame(
    Metabolite = c("Increased", "metab1", "Decreased", "metab11"),
    `Z-score` = c("", "2.45", "", "-1.51")
  )

  expect_silent(create_pdf_violin_plots(
    test_pdf_dir,
    test_patient_id,
    test_metab_perpage,
    test_top_metab_pt,
    test_explanation
  ))

  out_pdf_violinplots <- file.path(test_pdf_dir, "R_P2025M1.pdf")
  expect_true(file.exists(out_pdf_violinplots))
  content_pdf_violinplots <- pdftools::pdf_text(out_pdf_violinplots)
  expect_snapshot(content_pdf_violinplots)

  unlink(test_pdf_dir, recursive = TRUE)
})

testthat::test_that("create_violin_plot: Create a violin plot", {
  test_patient_id <- "P2025M1"
  test_sub_perpage <- "test acyl carnitines"

  test_acyl_carnitines_df <- read.delim(test_path("fixtures/", "test_acyl_carnitines_df.txt"))
  attr(test_acyl_carnitines_df, "y_order") <- rev(unique(test_acyl_carnitines_df$HMDB_name))

  test_patient_zscore_df <- test_acyl_carnitines_df %>% filter(Sample == test_patient_id)

  test_metab_zscores_df <- test_acyl_carnitines_df %>% filter(Sample != test_patient_id)

  expect_silent(create_violin_plot(test_metab_zscores_df, test_patient_zscore_df, test_sub_perpage, test_patient_id))

  expect_doppelganger("violin_plot_P2025M1", create_violin_plot(
    test_metab_zscores_df, test_patient_zscore_df,
    test_sub_perpage, test_patient_id
  ))
})

testthat::test_that("run_diem_algorithm :Run dIEM algorithm", {
  test_expected_biomarkers_df <- read.delim(test_path("fixtures/", "test_expected_biomarkers_df.txt"))
  test_zscore_patient_df <- read.delim(test_path("fixtures/", "test_zscore_patient_df.txt"))
  test_sample_cols <- c("P2025M1", "P2025M2", "P2025M3", "P2025M4")

  expect_equal(
    dim(run_diem_algorithm(test_expected_biomarkers_df, test_zscore_patient_df, test_sample_cols)),
    c(7, 5)
  )
  expect_identical(
    colnames(run_diem_algorithm(test_expected_biomarkers_df, test_zscore_patient_df, test_sample_cols)),
    c("Disease", "P2025M1", "P2025M2", "P2025M3", "P2025M4")
  )
  expect_identical(
    run_diem_algorithm(test_expected_biomarkers_df, test_zscore_patient_df, test_sample_cols)$Disease,
    c("Disease A", "Disease B", "Disease C", "Disease D", "Disease E", "Disease F", "Disease G")
  )
  expect_equal(run_diem_algorithm(test_expected_biomarkers_df, test_zscore_patient_df, test_sample_cols)$P2025M1,
    c(10.94172, 0.95343, 12.12121, 0.00000, 44.28850, 0.00000, -38.70370),
    tolerance = 0.0001
  )
})

testthat::test_that("rank_patient_zscores: Ranking Z-scores for a patient", {
  test_zscore_col <- c(1, 5, 6, 2, 7, -2, 3)

  expect_equal(length(rank_patient_zscores(test_zscore_col)), 7)

  expect_identical(
    rank_patient_zscores(test_zscore_col),
    c(6, 3, 2, 5, 1, 1, 4)
  )

  test_zscore_col <- c(3, 2, 1, 3)

  expect_identical(
    rank_patient_zscores(test_zscore_col),
    c(1, 2, 3, 1)
  )

  test_zscore_col <- c(-1, -2, -3, -4)

  expect_identical(
    rank_patient_zscores(test_zscore_col),
    c(4, 3, 2, 1)
  )
})

testthat::test_that("save_prob_scores_to_excel: Saving the probability score dataframe as an Excel file", {
  local_edition(3)
  test_probability_score_df <- read.delim(test_path("fixtures/", "test_probability_score_df.txt"))

  test_run_name <- "test_run"
  out_excel_file <- file.path(paste0("dIEM_algoritme_output_", test_run_name, ".xlsx"))

  expect_silent(save_prob_scores_to_excel(test_probability_score_df, test_run_name))
  expect_true(file.exists(out_excel_file))

  content_excel_file <- openxlsx::read.xlsx(out_excel_file, sheet = 1)
  expect_snapshot_output(content_excel_file)

  file.remove(out_excel_file, recursive = TRUE)
})

testthat::test_that("prepare_intensities_zscore_df: Preparing the intensities and Z-score dataframe", {
  test_intensities_zscore_df <- read.delim(test_path("fixtures/", "test_outlist.txt"))
  test_intensities_zscore_df$nr_ctrls <- c(28, 27, 30, 25)
  test_intensities_zscore_df$avg_ctrls <- c(16129.17, 1150, 1231.25, 4015.833)
  test_intensities_zscore_df$sd_ctrls <- c(51606.27, 349.6752, 313.7249, 6152.479)

  prepare_intensities_zscore_df(test_intensities_zscore_df)

  expect_equal(
    colnames(prepare_intensities_zscore_df(test_intensities_zscore_df)),
    c(
      "HMDB_code", "HMDB_name", "C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "C106.1", "C107.1", "C108.1",
      "C109.1", "C110.1", "C111.1", "C112.1", "P2.1", "P3.1", "mean_controls", "sd_controls"
    )
  )

  expect_equal(
    unname(sapply(prepare_intensities_zscore_df(test_intensities_zscore_df), class)),
    c(
      "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
      "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"
    )
  )
})

testthat::test_that("get_colnames_samples: Get all column names containing a specific prefix", {
  test_colnames <- c(
    "HMDB_name", "HMDB_code", "P1001", "P1001_Zscore", "P1002",
    "P1003", "P1004", "C101_Zscore", "C102", "C103"
  )
  test_intensities_zscore_df <- read.delim(test_path("fixtures/", "test_intensities_zscore_df.txt"))

  expect_equal(get_colnames_samples(test_intensities_zscore_df, "P"), c("P2025M1", "P2025M2", "P2025M3", "P2025M4", "P2025M5"))
  expect_equal(get_colnames_samples(test_intensities_zscore_df, "C"), c("C101.1", "C102.1", "C103.1", "C104.1", "C105.1"))
})

testthat::test_that("add_zscores_ratios_to_df: Add Zscores for multiple ratios to the dataframe", {
  test_outlist_df <- read.delim(test_path("fixtures/GenerateViolinPlots", "test_outlist_df.txt"))

  test_metabolites_ratios_df <- data.frame(
    HMDB.code = c("HMDB000TT1", "HMDB000TT2", "HMDB000TT3"),
    Ratio_name = c("Test_ratio1", "Test_ratio2", "Test_ratio3"),
    HMDB_numerator = c("HMDB001", "HMDB003plusHMDB003", "HMDB001plusHMDB003plusHMDB004"),
    HMDB_denominator = c("HMDB002", "HMDB004", "one")
  )

  test_all_sample_ids <- c(
    "C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "P2025M1", "P2025M2", "P2025M3",
    "P2025M4", "P2025M5"
  )

  expect_equal(
    add_zscores_ratios_to_df(test_outlist_df, test_metabolites_ratios_df, test_all_sample_ids)$HMDB_code,
    c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3")
  )
  expect_equal(
    add_zscores_ratios_to_df(
      test_outlist_df,
      test_metabolites_ratios_df,
      test_all_sample_ids
    )$C101.1,
    c(1000, 1200, 1300, 1400, 1500, 1600, -0.2630344, -0.1069152, 11.8533096)
  )
  expect_equal(
    add_zscores_ratios_to_df(
      test_outlist_df,
      test_metabolites_ratios_df,
      test_all_sample_ids
    )$C101.1_Zscore,
    c(0.45, 1.67, -1.86, 0.58, 2.47, -0.56, -0.5899371, 0.4858991, -0.4552026),
    tolerance = 0.0001
  )
})

testthat::test_that("calculate_zscore_ratios: Calculate Zscores for ratios", {
  test_outlist_df <- read.delim(test_path("fixtures/GenerateViolinPlots", "test_outlist_df.txt"))
  test_intensities_zscores_df <- prepare_intensities_zscore_df(test_outlist_df)

  test_metabolites_ratios_df <- data.frame(
    HMDB.code = c("HMDB000TT1", "HMDB000TT2", "HMDB000TT3"),
    Ratio_name = c("Test_ratio1", "Test_ratio2", "Test_ratio3"),
    HMDB_numerator = c("HMDB001", "HMDB003plusHMDB003", "HMDB001plusHMDB003plusHMDB004"),
    HMDB_denominator = c("HMDB002", "HMDB004", "one")
  )

  test_all_sample_ids <- c(
    "C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "P2025M1", "P2025M2", "P2025M3",
    "P2025M4", "P2025M5"
  )

  expect_equal(
    calculate_zscore_ratios(test_metabolites_ratios_df, test_outlist_df, test_all_sample_ids)$HMDB_code,
    c("HMDB000TT1", "HMDB000TT2", "HMDB000TT3")
  )
  expect_equal(
    calculate_zscore_ratios(test_metabolites_ratios_df, test_outlist_df, test_all_sample_ids)$C101.1,
    c(-0.2630344, -0.1069152, 11.8533096)
  )
  expect_equal(
    calculate_zscore_ratios(test_metabolites_ratios_df, test_outlist_df, test_all_sample_ids)$C101.1_Zscore,
    c(-0.5899371, 0.4858991, -0.4552026),
    tolerance = 0.0001
  )
})

testthat::test_that("get_list_page_plot_data: Get a list of dataframes for each chunk", {
  test_metabolite_class_patients_df <- read.delim(test_path("fixtures/GenerateViolinPlots",
                                                            "test_metabolite_class_patients_df.txt"))
  test_metabolite_class_controls_df <- read.delim(test_path("fixtures/GenerateViolinPlots",
                                                            "test_metabolite_class_controls_df.txt"))
  test_nr_plots_perpage <- 2
  test_metabolite_in_chunks <- list(
    c("HMDB001", "HMDB002", "HMDB003"),
    c("HMDB004", "HMDB011", "HMDB012"),
    c("HMDB000TT1", "HMDB000TT2", "HMDB000TT3")
  )
  
  t <- get_list_page_plot_data(
    test_metabolite_in_chunks,
    test_metabolite_class_patients_df,
    test_metabolite_class_controls_df,
    test_nr_plots_perpage
  )
  
  expect_type(get_list_page_plot_data(
    test_metabolite_in_chunks,
    test_metabolite_class_patients_df,
    test_metabolite_class_controls_df,
    test_nr_plots_perpage
  ), "list")
  
  expect_equal(length(get_list_page_plot_data(
    test_metabolite_in_chunks,
    test_metabolite_class_patients_df,
    test_metabolite_class_controls_df,
    test_nr_plots_perpage
  )), 3)
  
  test_plot_data_list <- get_list_page_plot_data(
    test_metabolite_in_chunks,
    test_metabolite_class_patients_df,
    test_metabolite_class_controls_df,
    test_nr_plots_perpage
  )
  for (num_chunk in seq_along(test_metabolite_in_chunks)) {
    expect_equal(unique(test_plot_data_list[[num_chunk]]$HMDB_name), test_metabolite_in_chunks[[num_chunk]])
    expect_equal(unique(test_plot_data_list[[num_chunk]]$Sample),
                 c("P2025M1", "P2025M2", "P2025M3", "C101.1", "C102.1", "C103.1"))
  }
  
})

testthat::test_that("make_and_save_violin_plot_pdfs: Make and save violin plots for each patient in a PDF", {
  test_zscore_patients_df <- read.delim(test_path("fixtures/GenerateViolinPlots", "test_zscore_patients_df.txt"))
  test_zscore_controls_df <- read.delim(test_path("fixtures/GenerateViolinPlots", "test_zscore_controls_df.txt"))
  test_path_metabolite_groups <- test_path("fixtures/test_metabolite_groups")
  test_nr_plots_perpage <- 2
  test_number_of_samples <- list(
    controls = 5,
    patients = 5
  )
  test_run_name <- "unit_test"
  test_protocol_name <- "UNIT_TEST_PROTOCOL"
  test_explanation_violin_plot <- "Unit test violin plot pdfs"
  test_number_of_metabolites <- list(
    highest = 2,
    lowest = 1
  )

  expect_silent(make_and_save_violin_plot_pdfs(
    test_zscore_patients_df,
    test_zscore_controls_df,
    test_path_metabolite_groups,
    test_nr_plots_perpage,
    test_number_of_samples,
    test_run_name,
    test_protocol_name,
    test_explanation_violin_plot,
    test_number_of_metabolites
  ))
  
  patient_ids <- c("P2025M1", "P2025M2", "P2025M3", "P2025M4", "P2025M5")
  for (patient_id in patient_ids) {
    pdf_file_name_diagnotics <- file.path(paste0("Diagnostics/MB", gsub("^P|M", "", patient_id), "_DIMS_PL_DIAG.pdf"))
    pdf_file_name_other <- file.path(paste0("Other/R_", patient_id, ".pdf"))
    
    expect_true(file.exists(pdf_file_name_diagnotics))
    expect_true(file.exists(pdf_file_name_other))
  }
  expect_true(file.exists("output_Helix_unit_test.csv"))
  
  unlink(c("Diagnostics", 'Other'), recursive = TRUE)
  file.remove("output_Helix_unit_test.csv")
})

testthat::test_that("get_probabilities_top_iems: Get the IEM probabilities for a patient for all diseases", {
  test_patient_top_iems_probs <- data.frame(
    Disease = c("Disease A", "Disease B", "Disease C", "Disease D"),
    P2025M1 = c(100, 75, 50, 25)
  )
  test_expected_biomarkers_df <- read.delim(test_path("fixtures", "test_expected_biomarkers_df.txt"))
  test_patient_id <- "P2025M1"
  

  expect_type(get_probabilities_top_iems(test_patient_top_iems_probs, test_expected_biomarkers_df, test_patient_id),
              "list")
  
  expect_equal(length(get_probabilities_top_iems(test_patient_top_iems_probs, test_expected_biomarkers_df, test_patient_id)),
               4)
  expect_equal(names(get_probabilities_top_iems(test_patient_top_iems_probs, test_expected_biomarkers_df, test_patient_id)),
               c("Disease A, probability score 100", "Disease B, probability score 75", "Disease C, probability score 50",
                 "Disease D, probability score 25"))
  
  expect_equal(get_probabilities_top_iems(test_patient_top_iems_probs,
                                          test_expected_biomarkers_df,
                                          test_patient_id)$"Disease A, probability score 100",
               data.frame(
                 HMDB_code = c("HMDB002", "HMDB012"),
                 HMDB_name = c("metab2", "metab12")
               ))
  
})

testthat::test_that("make_and_save_diem_plots: Make and save dIEM plots", {
  test_diem_probability_score <- data.frame(
    Disease = c("Disease A", "Disease B", "Disease C", "Disease D", "Disease E"),
    P2025M1 = c(100, 75, 50, 25, 12.5),
    P2025M2 = c(25, 0, 2, 8, 3),
    P2025M3 = c(0, 1, 2, 3, 4)
  )
  test_patient_ids <- c("P2025M1", "P2025M2", "P2025M3")
  test_expected_biomarkers_df <- read.delim(test_path("fixtures", "test_expected_biomarkers_df.txt"))
  test_zscore_patients_df <- read.delim(test_path("fixtures/GenerateViolinPlots", "test_zscore_patients_df.txt"))
  test_zscore_controls_df <- read.delim(test_path("fixtures/GenerateViolinPlots", "test_zscore_controls_df.txt"))
  test_nr_plots_perpage <- 2
  test_number_of_samples <- list(
    controls = 5,
    patients = 5
  )
  test_number_of_metabolites <- list(
    highest = 2,
    lowest = 1
  )
  test_iem_variables <- list(
    top_number_iem_diseases = 5,
    threshold_iem = 5
  )
  test_explanation_violin_plot <- "Unit test violin plot pdfs"

  expect_silent(make_and_save_diem_plots(
    test_diem_probability_score,
    test_patient_ids,
    test_expected_biomarkers_df,
    test_zscore_patients_df,
    test_zscore_controls_df,
    test_nr_plots_perpage,
    test_number_of_samples,
    test_number_of_metabolites,
    test_iem_variables,
    test_explanation_violin_plot
  ))
  unlink("dIEM_plots/", recursive = TRUE)
  
  expect_equal(make_and_save_diem_plots(
    test_diem_probability_score,
    test_patient_ids,
    test_expected_biomarkers_df,
    test_zscore_patients_df,
    test_zscore_controls_df,
    test_nr_plots_perpage,
    test_number_of_samples,
    test_number_of_metabolites,
    test_iem_variables,
    test_explanation_violin_plot
  ), "P2025M3")
  
  expect_true(file.exists("dIEM_plots/IEM_P2025M1.pdf"))
  expect_true(file.exists("dIEM_plots/IEM_P2025M2.pdf"))
  
  unlink("dIEM_plots/", recursive = TRUE)
})

testthat::test_that("make_metabolite_order: Make the order of metabolites for the violin plots", {
  test_metabolites_vector <- c("metab1", "metab2", "metab3")
  test_num_plots_per_page <- 5
  
  make_metabolite_order(test_num_plots_per_page, test_metabolites_vector)
  
  expect_equal(length(make_metabolite_order(test_num_plots_per_page, test_metabolites_vector)), 5)
  expect_equal(make_metabolite_order(test_num_plots_per_page, test_metabolites_vector),
               c("metab1", "metab2", "metab3", "  ", "   "))
})

testthat::test_that("pad_truncate_hmdb_names: Pad or truncate HMDB names to a fixed width", {
  test_dataframe <- data.frame(
    HMDB_name = c("metab1", "metabolite2", "metabo3")
  )
  test_width <- 10
  test_pad_character <- "+"
  
  expect_equal(nrow(pad_truncate_hmdb_names(test_dataframe, test_width, test_pad_character)),
               3)
  expect_equal(pad_truncate_hmdb_names(test_dataframe, test_width, test_pad_character)$HMDB_name,
               c("metab1++++", "metabol...", "metabo3+++"))
  
  test_df <- pad_truncate_hmdb_names(test_dataframe, test_width, test_pad_character)
  for (row in seq(nrow(test_df))) {
    expect_equal(nchar(test_df[row, "HMDB_name"]), 10)
  }
})

testthat::test_that("save_patient_no_iem: Save a list of patient IDs to a text file", {
  local_edition(3)
  test_threshold_iem <- 5
  test_patient_no_iem <- c("Patient1", "Patient2")
  
  expect_silent(save_patient_no_iem(test_threshold_iem, test_patient_no_iem))
  
  expect_true(file.exists("missing_probability_scores.txt"))
  
  expect_snapshot_file("missing_probability_scores.txt")
  file.remove("missing_probability_scores.txt")
})
