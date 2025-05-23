# unit tests for CheckQC
# functions: get_internal_standards, save_internal_standard_plot,
#            get_pos_ctrl_data, round_df,
#            get_is_intensities, calculate_coefficient_of_variation,
#            check_missing_mz
library(ggplot2)
suppressMessages(library("dplyr"))
source("../../export/check_qc_functions.R")

testthat::test_that("Get the internal standards data", {
  test_internal_standards <- read.delim(test_path("fixtures", "test_internal_standards.txt"))
  test_outlist <- read.delim(test_path("fixtures", "test_outlist_IS.txt"))

  repl_pattern <- list(
    C101.1 = c("test1", "test2"),
    C102.1 = c("test3", "test4"),
    P2.1 = c("test5", "test6"),
    P3.1 = c("test7", "test8")
  )

  expect_equal(colnames(get_internal_standards(test_internal_standards, "summed", repl_pattern, "Plasma", "date_test", "Test")),
               c("HMDB_code", "HMDB_name", "Sample", "Intensity", "Matrix", "Rundate", "Project"))
  expect_equal(nrow(get_internal_standards(test_internal_standards, "summed", repl_pattern, "Plasma", "date_test", "Test")), 16)
  expect_equal(get_internal_standards(test_internal_standards, "summed", repl_pattern, "Plasma", "date_test", "Test")$Intensity,
               c(100, 200, 300, 400, 125, 225, 325, 425, 150, 250, 350, 450, 175, 275, 375, 475))

  expect_equal(colnames(get_internal_standards(test_internal_standards, "pos", test_outlist, "Plasma", "date_test", "Test")),
               c("HMDB_code", "HMDB_name", "Sample", "Intensity", "Matrix", "Rundate", "Project"))
  expect_equal(nrow(get_internal_standards(test_internal_standards, "pos", test_outlist, "Plasma", "date_test", "Test")), 8)
  expect_equal(get_internal_standards(test_internal_standards, "pos", test_outlist, "Plasma", "date_test", "Test")$Intensity,
               c(100, 300, 125, 325, 150, 350, 175, 375))
  expect_equal(unique(get_internal_standards(test_internal_standards, "pos", test_outlist, "Plasma", "date_test", "Test")$HMDB_code),
               c("HMDB100000", "HMDB300000"))
})

testthat::test_that("Save internal standard plots", {
  local_edition(3)
  temp_dir <- "./"
  dir.create(paste0(temp_dir, "plots"))

  test_plot_data <- data.frame(
    HMDB_code = rep(c("HMDB100000", "HMDB200000"), each = 3),
    HMDB_name = rep(c("metab_1 (IS)", "metab_2 (IS)"), each = 3),
    Sample = c("C101.1", "C102.1", "P2.1", "C101.1", "C102.1", "P2.1"),
    Intensity = c(100, 300, 125, 325, 150, 350),
    Matrix = rep("Plasma", each = 6),
    Rundate = rep("Date_test", each = 6),
    Project = rep("Test", each = 6)
  )
  test_plot_data$Sample_level <- factor(test_plot_data$Sample, levels = unique(test_plot_data$Sample))

  test_hline_data <- data.frame(
    int_line = c(105, 250),
    HMDB_name = c("metab_1 (IS)", "metab_2 (IS)")
  )

  file_name_barplot <- "test_barplot"
  out_file_barplot <- file.path(temp_dir, "/plots", paste0(file_name_barplot, ".png"))

  expect_silent(
    save_internal_standard_plot(
      test_plot_data, "barplot", "Test barplot", temp_dir,
      file_name_barplot, 6, 4, test_hline_data
    )
  )
  expect_true(file.exists(out_file_barplot))
  expect_snapshot_file(out_file_barplot, "test_barplot.png")

  file_name_lineplot <- "test_lineplot"
  out_file_lineplot <- file.path(temp_dir, "plots", paste0(file_name_lineplot, ".png"))
  expect_silent(
    save_internal_standard_plot(
      test_plot_data, "lineplot", "Test lineplot", temp_dir,
      file_name_lineplot, 6, 4
    )
  )
  expect_true(file.exists(out_file_lineplot))

  test_plot_data <- test_plot_data[-c(1, 2, 3, 4, 5, 6), ]

  expect_identical(save_internal_standard_plot(
    test_plot_data, "barplot", "Test barplot", temp_dir,
    file_name_barplot, 6, 4, test_hline_data
  ), NULL)

  file_name_barplot_select <- "select_test_barplot"

  expect_silent(
    save_internal_standard_plot(
      test_plot_data, "barplot", "Test barplot", temp_dir,
      file_name_barplot_select, 6, 4, test_hline_data
    )
  )
  unlink("plots", recursive = TRUE)
})

testthat::test_that("Get positive control data", {
  test_outlist <- data.frame(
    plots = c(NA, NA, NA, NA),
    C101.1 = c(100, 200, 300, 400),
    C102.1 = c(125, 225, 325, 425),
    P1002.1 = c(150, 250, 350, 450),
    P1003.1 = c(175, 275, 375, 475),
    HMDB_name = c("Propionylcarnitine", "Propionylglycine", "Glycine", "L-Phenylalanine"),
    HMDB_key = c("HMDB0000824", "HMDB0000725", "HMDB0000123", "HMDB0000159"),
    HMDB_code = c("HMDB0000824", "HMDB0000725", "HMDB0000123", "HMDB0000159"),
    name = c("Propionylcarnitine", "Propionylglycine", "Glycine", "L-Phenylalanine"),
    avg.ctrls = c(137.5, 237.5, 337.5, 437.5),
    sd.ctrls = c(32.2749, 32.2749, 32.2749, 32.2749),
    nr.ctrls = c(2, 2, 2, 2),
    C101.1_Zscore = c(0.25, 1.65, 2.75, -0.35),
    C102.2_Zscore = c(1.89, -0.42, 0.22, 1.11),
    P1002.1_Zscore = c(0.59, 1.36, -0.51, -0.32),
    P1003.1_Zscore = c(1.03, 0.28, 0.78, 0.68)
  )
  rownames(test_outlist) <- c("HMDB0000824", "HMDB0000725", "HMDB0000123", "HMDB0000159")

  pa_sample_name_test <- "P1002.1_Zscore"
  pa_codes_test <- c("HMDB0000824", "HMDB0000725", "HMDB0000123")
  pa_names_test <- c("Propionylcarnitine", "Propionylglycine", "Glycine")

  expect_identical(colnames(get_pos_ctrl_data(test_outlist, pa_sample_name_test, pa_codes_test, pa_names_test)),
                   c("HMDB_code", "HMDB_name", "Sample", "Zscore"))
  expect_equal(get_pos_ctrl_data(test_outlist, pa_sample_name_test, pa_codes_test, pa_names_test)$Zscore,
               c(0.59, 1.36, -0.51))
  expect_identical(get_pos_ctrl_data(test_outlist, pa_sample_name_test, pa_codes_test, pa_names_test)$HMDB_code,
                   c("HMDB0000824", "HMDB0000783", "HMDB0000123"))

  pku_sample_name_test <- "P1003.1_Zscore"
  pku_codes_test <- c("HMDB0000159")
  pku_names_test <- c("L-Phenylalanine")

  expect_identical(colnames(get_pos_ctrl_data(test_outlist, pku_sample_name_test, pku_codes_test, pku_names_test)),
                   c("HMDB_code", "HMDB_name", "Sample", "Zscore"))
  expect_equal(get_pos_ctrl_data(test_outlist, pku_sample_name_test, pku_codes_test, pku_names_test)$Zscore,
               0.68)
  expect_identical(get_pos_ctrl_data(test_outlist, pku_sample_name_test, pku_codes_test, pku_names_test)$HMDB_code,
                   "HMDB0000159")
})

testthat::test_that("Rounding numbers in a dataframe", {
  test_df <- data.frame(
    C101.1 = c(1.54641, 2.457, -0.4568, 0.6468),
    C102.1 = c(-0.4685, 0.5986, 2.836, -1.156),
    HMDB_name = c("metab_1 (IS)", "metab_85", "metab_3 (IS)", "metab_245"),
    HMDB_ID_all = c("HMDB100000", "HMDB68425", "HMDB300000", "HMDB84684")
  )

  expect_identical(colnames(round_df(test_df, 0)), c("C101.1", "C102.1", "HMDB_name", "HMDB_ID_all"))
  expect_equal(round_df(test_df, 2)$C101.1, c(1.55, 2.46, -0.46, 0.65), tolerance = 0.001)
  expect_equal(round_df(test_df, 4)$C101.1, c(1.5464, 2.4570, -0.4568, 0.6468), tolerance = 0.001)
  expect_equal(round_df(test_df, 0)$C101.1, c(2, 2, 0, 1), tolerance = 0.001)
})

testthat::test_that("Calculate coefficient of variation", {
  test_internal_standards <- data.frame(
    C101.1 = c(100.5, 200.54, 300.15, 400.1),
    C102.1 = c(125, 225, 325, 425),
    P2.1 = c(150, 250, 350, 450),
    P3.1 = c(175, 275, 375, 475)
  )
  rownames(test_internal_standards) <- c("HMDB100000", "HMDB200000", "HMDB300000", "HMDB400000")

  expect_identical(colnames(calculate_coefficient_of_variation(test_internal_standards)),
                   c("CV_perc", "mean", "sd", "C101.1", "C102.1", "P2.1", "P3.1"))
  expect_equal(calculate_coefficient_of_variation(test_internal_standards)$CV_perc,
               c(23.2, 13.4, 9.5, 7.3))
  expect_equal(calculate_coefficient_of_variation(test_internal_standards)$mean,
               c(138, 238, 338, 438))
})

testthat::test_that("Get internal standard intensities", {
  test_outlist_IS <- read.delim(test_path("fixtures", "test_outlist_IS.txt"))
  test_internal_standards <- test_outlist_IS[grepl("IS", test_outlist_IS$HMDB_name), ]
  
  int_cols <- c(1, 2, 3, 4)

  expect_identical(colnames(get_is_intensities(test_internal_standards, int_cols = int_cols)),
                   c("IS_name", "CV_perc", "mean", "sd", "C101.1", "C102.1", "P2.1", "P3.1"))
  expect_equal(get_is_intensities(test_internal_standards, int_cols = int_cols)$CV_perc,
               c(23.2, 9.5))

  is_codes <- c("HMDB100000", "HMDB300000")

  expect_identical(rownames(get_is_intensities(test_internal_standards, is_codes = is_codes)),
                   c("HMDB100000", "HMDB300000"))
  expect_equal(get_is_intensities(test_internal_standards, is_codes = is_codes)$C101.1,
               c(100, 300))
})

testthat::test_that("Check missing mz values", {
  test_mz_pgrp <- seq(70, 599, by = 1)

  expect_identical(check_missing_mz(test_mz_pgrp, "Test"),
                   "Test mode did not have missing mz values")

  test_mz_pgrp <- test_mz_pgrp[! test_mz_pgrp %in% seq(550, 560, by = 1)]

  expect_identical(check_missing_mz(test_mz_pgrp, "Test")[[1]],
                   "Missing m/z values Test mode")
  expect_identical(check_missing_mz(test_mz_pgrp, "Test")$`1`,
                   c(550, 551, 552, 553, 554, 555, 556, 557, 558, 559, 560))
})
