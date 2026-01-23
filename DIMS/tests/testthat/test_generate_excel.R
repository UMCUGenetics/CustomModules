# unit tests for GenerateExcel
# functions: get_intensities_cols, calculate_zscores,
#            robust_scaler, remove_outliers_grubbs,
#            save_to_rdata_and_txt, set_row_height_col_width_wb
library("ggplot2")
library("reshape2")
library("openxlsx")
library("vdiffr")
suppressMessages(library("tidyr"))
suppressMessages(library("dplyr"))
suppressMessages(library("stringr"))
source("../../export/generate_excel_functions.R")

testthat::test_that("get_intensities_cols: Get indices of columns and dataframe of intensities of a given label", {
  test_outlist <- read.delim(test_path("fixtures", "test_outlist.txt"))

  control_label <- "C"
  case_label <- "P"

  expect_equal(get_intensities_cols(test_outlist, control_label)$col_idx, c(2:13))
  expect_equal(
    colnames(get_intensities_cols(test_outlist, control_label)$df_intensities),
    c("C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "C106.1", "C107.1", "C108.1", "C109.1", "C110.1", "C111.1", "C112.1")
  )
  expect_equal(
    rownames(get_intensities_cols(test_outlist, control_label)$df_intensities),
    c("HMDB001", "HMDB002", "HMDB003", "HMDB004")
  )
  expect_equal(get_intensities_cols(test_outlist, control_label)$df_intensities$C101.1, c(1000, 1200, 1300, 1400))

  expect_equal(get_intensities_cols(test_outlist, case_label)$col_idx, c(14, 15))
  expect_equal(colnames(get_intensities_cols(test_outlist, case_label)$df_intensities), c("P2.1", "P3.1"))
})

testthat::test_that("calculate_zscores: Calculating Z-scores using different methods for excluding controls", {
  test_outlist <- read.delim(test_path("fixtures", "test_outlist.txt"))
  control_intensities <- read.delim(test_path("fixtures", "test_control_intensities.txt"))

  control_col_idx <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
  intensity_col_ids <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
  startcol <- ncol(test_outlist) + 4
  perc <- 5
  outlier_threshold <- 2

  expect_type(
    calculate_zscores(
      test_outlist,
      "_Zscore",
      control_intensities,
      NULL,
      intensity_col_ids,
      startcol
    ),
    "list"
  )

  expect_identical(
    colnames(
      calculate_zscores(
        test_outlist,
        "_Zscore",
        control_intensities,
        NULL,
        intensity_col_ids,
        startcol
      )
    ),
    c(
      "plots", "C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "C106.1", "C107.1", "C108.1", "C109.1", "C110.1", "C111.1",
      "C112.1", "P2.1", "P3.1", "HMDB_name", "HMDB_name_all", "HMDB_ID_all", "sec_HMDB_ID", "HMDB_key", "sec_HMDB_ID_rlvc",
      "name", "relevance", "descr", "origin", "fluids", "tissue", "disease", "pathway", "HMDB_code", "avg_ctrls", "sd_ctrls",
      "nr_ctrls", "C101.1_Zscore", "C102.1_Zscore", "C103.1_Zscore", "C104.1_Zscore", "C105.1_Zscore", "C106.1_Zscore",
      "C107.1_Zscore", "C108.1_Zscore", "C109.1_Zscore", "C110.1_Zscore", "C111.1_Zscore", "C112.1_Zscore", "P2.1_Zscore",
      "P3.1_Zscore"
    )
  )
  expect_equal(
    round(
      calculate_zscores(
        test_outlist,
        "_Zscore",
        control_intensities,
        NULL,
        intensity_col_ids,
        startcol
      )$avg_ctrls, 3
    ),
    c(16129.167, 1150.0, 1231.250, 4015.833),
    tolerance = 0.001
  )
  expect_equal(
    calculate_zscores(
      test_outlist,
      "_Zscore",
      control_intensities,
      NULL,
      intensity_col_ids,
      startcol
    )$P2.1_Zscore,
    c(-0.2544103, 32.4586955, 13.6066674, 0.4037668),
    tolerance = 0.001
  )

  expect_type(
    calculate_zscores(
      test_outlist,
      "_RobustZscore",
      control_col_idx,
      perc,
      intensity_col_ids,
      startcol
    ),
    "list"
  )

  expect_identical(
    colnames(
      calculate_zscores(
        test_outlist,
        "_RobustZscore",
        control_col_idx,
        perc,
        intensity_col_ids,
        startcol
      )
    )[34:47],
    c(
      "C101.1_RobustZscore", "C102.1_RobustZscore", "C103.1_RobustZscore", "C104.1_RobustZscore", "C105.1_RobustZscore",
      "C106.1_RobustZscore", "C107.1_RobustZscore", "C108.1_RobustZscore", "C109.1_RobustZscore", "C110.1_RobustZscore",
      "C111.1_RobustZscore", "C112.1_RobustZscore", "P2.1_RobustZscore", "P3.1_RobustZscore"
    )
  )

  expect_equal(
    calculate_zscores(
      test_outlist,
      "_RobustZscore",
      control_col_idx,
      perc,
      intensity_col_ids,
      startcol
    )$avg_ctrls,
    c(1255.0, 1110.0, 1227.5, 2811.5),
    tolerance = 0.001
  )

  expect_equal(
    calculate_zscores(
      test_outlist,
      "_RobustZscore",
      control_col_idx,
      perc,
      intensity_col_ids,
      startcol
    )$P2.1_RobustZscore,
    c(9.1511750, 46.9804468, 16.8039663, 0.8565111),
    tolerance = 0.001
  )

  expect_type(
    calculate_zscores(
      test_outlist,
      "_OutlierRemovedZscore",
      control_col_idx,
      outlier_threshold,
      intensity_col_ids,
      startcol
    ),
    "list"
  )

  expect_identical(
    colnames(
      calculate_zscores(
        test_outlist,
        "_OutlierRemovedZscore",
        control_col_idx,
        outlier_threshold,
        intensity_col_ids,
        startcol
      )
    )[34:47],
    c(
      "C101.1_OutlierRemovedZscore", "C102.1_OutlierRemovedZscore", "C103.1_OutlierRemovedZscore",
      "C104.1_OutlierRemovedZscore", "C105.1_OutlierRemovedZscore", "C106.1_OutlierRemovedZscore",
      "C107.1_OutlierRemovedZscore", "C108.1_OutlierRemovedZscore", "C109.1_OutlierRemovedZscore",
      "C110.1_OutlierRemovedZscore", "C111.1_OutlierRemovedZscore", "C112.1_OutlierRemovedZscore",
      "P2.1_OutlierRemovedZscore", "P3.1_OutlierRemovedZscore"
    )
  )

  expect_equal(
    calculate_zscores(
      test_outlist,
      "_OutlierRemovedZscore",
      control_col_idx,
      outlier_threshold,
      intensity_col_ids,
      startcol
    )$avg_ctrls,
    c(1231.818, 1077.273, 1231.250, 2649.091),
    tolerance = 0.001
  )
  expect_equal(
    calculate_zscores(
      test_outlist,
      "_OutlierRemovedZscore",
      control_col_idx,
      outlier_threshold,
      intensity_col_ids,
      startcol
    )$nr_ctrls,
    c(11, 11, 12, 11)
  )
  expect_equal(
    calculate_zscores(
      test_outlist,
      "_OutlierRemovedZscore",
      control_col_idx,
      outlier_threshold,
      intensity_col_ids,
      startcol
    )$P2.1_OutlierRemovedZscore,
    c(8.9955723, 44.9136860, 13.6066674, 0.9345077),
    tolerance = 0.001
  )
})

testthat::test_that("robust_scaler: Use robust scaler", {
  control_intensities <- read.delim(test_path("fixtures", "test_control_intensities.txt"))

  control_col_idx <- c(1)
  perc <- 5

  expect_type(robust_scaler(control_intensities[1, 1:12], control_col_idx, perc), "double")
  expect_length(robust_scaler(control_intensities[1, 1:12], control_col_idx, perc), 10)
  expect_equal(
    robust_scaler(
      control_intensities[1, 1:12],
      control_col_idx,
      perc
    ),
    c(1050, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1450, 1650),
    tolerance = 0.001
  )
})

testthat::test_that("remove_outliers_grubbs: Use Grubbs outlier removal", {
  control_intensities <- read.delim(test_path("fixtures", "test_control_intensities.txt"))

  outlier_threshold <- 2

  expect_type(remove_outliers_grubbs(control_intensities[1, ], outlier_threshold), "list")
  expect_length(remove_outliers_grubbs(control_intensities[1, ], outlier_threshold), 11)
  expect_equal(
    as.numeric(remove_outliers_grubbs(
      control_intensities[1, ],
      outlier_threshold
    )),
    c(1000, 1100, 1300, 1650, 1050, 1150, 1350, 1450, 1200, 1050, 1250),
    tolerance = 0.001
  )
  expect_identical(
    colnames(
      remove_outliers_grubbs(
        control_intensities[1, ],
        outlier_threshold
      )
    ),
    c("C101.1", "C102.1", "C103.1", "C104.1", "C106.1", "C107.1", "C108.1", "C109.1", "C110.1", "C111.1", "C112.1")
  )
})

testthat::test_that("save_to_rdata_and_txt: Save data to RData and txt file", {
  test_df <- data.frame(
    C101.1 = c(100, 200, 300, 400),
    C102.1 = c(125, 225, 325, 425),
    P2.1 = c(150, 250, 350, 450),
    P3.1 = c(175, 275, 375, 475),
    HMDB_name = c("metab_1", "metab_2", "metab_3", "metab_4"),
    sec_HMDB_ID = c("HMDB1;HMDB11", "", "HMDB3;HMDB13", "HMDB4")
  )
  rownames(test_df) <- c("HMDB001", "HMDB002", "HMDB003", "HMDB004")

  expect_no_error(save_to_rdata_and_txt(test_df, "test_df"))
  expect_true(all(file.exists(c("test_df.txt", "test_df.RData"))))

  check_cols_and_values <- function(df) {
    expect_identical(colnames(df), c("C101.1", "C102.1", "P2.1", "P3.1", "HMDB_name", "sec_HMDB_ID"))
    expect_equal(df$P2.1, c(150, 250, 350, 450))
  }

  check_cols_and_values(read.delim("test_df.txt"))
  check_cols_and_values(get(load("test_df.RData")))

  file.remove("test_df.txt", "test_df.RData")
})

testthat::test_that("set_row_height_col_width_wb: Check row height and column width in a workbook", {
  test_wb_plots <- openxlsx::createWorkbook("Test")
  openxlsx::addWorksheet(test_wb_plots, "Test_with_plots")

  sheetname_with_plots <- "Test_with_plots"
  num_rows_df <- 5
  num_cols_df <- 5
  plot_width <- 100

  correct_col_widths <- c("5", "20", "20", "20", "20")
  names(correct_col_widths) <- c(1, 2, 3, 4, 5)
  attr(correct_col_widths, "hidden") <- c("0", "0", "0", "0", "0")
  expect_identical(
    set_row_height_col_width_wb(
      test_wb_plots,
      sheetname_with_plots,
      num_rows_df,
      num_cols_df,
      plot_width,
      plots_present = TRUE
    )$colWidths[[1]],
    correct_col_widths
  )

  correct_row_heights <- c("140", "140", "140", "140", "140")
  names(correct_row_heights) <- c(2, 3, 4, 5, 6)
  expect_identical(
    set_row_height_col_width_wb(test_wb_plots,
      sheetname_with_plots,
      num_rows_df,
      num_cols_df,
      plot_width,
      plots_present = TRUE
    )$rowHeights[[1]],
    correct_row_heights
  )

  rm(test_wb_plots)

  test_wb_no_plots <- openxlsx::createWorkbook("Test")
  openxlsx::addWorksheet(test_wb_no_plots, "Test_no_plots")
  sheetname_no_plots <- "Test_no_plots"

  correct_col_widths <- c("20", "20", "20", "20", "20")
  names(correct_col_widths) <- c(1, 2, 3, 4, 5)
  attr(correct_col_widths, "hidden") <- c("0", "0", "0", "0", "0")
  expect_identical(
    set_row_height_col_width_wb(
      test_wb_no_plots,
      sheetname_no_plots,
      num_rows_df,
      num_cols_df,
      plot_width = NULL,
      plots_present = FALSE
    )$colWidths[[1]],
    correct_col_widths
  )

  correct_row_heights <- c("18", "18", "18", "18", "18")
  names(correct_row_heights) <- c(1, 2, 3, 4, 5)
  expect_identical(
    set_row_height_col_width_wb(
      test_wb_no_plots,
      sheetname_no_plots,
      num_rows_df,
      num_cols_df,
      plot_width = NULL,
      plots_present = FALSE
    )$rowHeights[[1]],
    correct_row_heights
  )

  rm(test_wb_no_plots)
})

testthat::test_that("transform_ints_df_plots: Check transformation of dataframe to long format", {
  test_intensities_plots_df <- data.frame(
    C101.1 = c(100, 200, 300, 400),
    C102.1 = c(125, 225, 325, 425),
    P2000M00002.1 = c(150, 250, 350, 450),
    P3000M00003.1 = c(175, 275, 375, 475),
    HMDB_key = c("metab_1", "metab_2", "metab_3", "metab_4")
  )

  test_row_index <- 1

  expect_equal(dim(transform_ints_df_plots(test_intensities_plots_df, test_row_index)), c(4, 4))
  expect_identical(
    colnames(
      transform_ints_df_plots(
        test_intensities_plots_df,
        test_row_index
      )
    ),
    c("Samples", "Intensities", "type", "group_size")
  )
  expect_identical(
    transform_ints_df_plots(
      test_intensities_plots_df,
      test_row_index
    )$Samples,
    c("C", "C", "P2000M00002", "P3000M00003")
  )
  expect_identical(
    transform_ints_df_plots(
      test_intensities_plots_df,
      test_row_index
    )$Intensities,
    c(100, 125, 150, 175)
  )
  expect_identical(
    transform_ints_df_plots(
      test_intensities_plots_df,
      test_row_index
    )$type,
    c("Control", "Control", "Patients", "Patients")
  )
  expect_equal(
    transform_ints_df_plots(
      test_intensities_plots_df,
      test_row_index
    )$group_size,
    c(2, 2, 1, 1)
  )

  test_row_index <- 2
  expect_identical(
    transform_ints_df_plots(
      test_intensities_plots_df,
      test_row_index
    )$Intensities,
    c(200, 225, 250, 275)
  )
})

testthat::test_that("create_boxplot: Create a boxplot for the Excel", {
  test_ints_df_long <- data.frame(
    Samples = c("C", "C", "C", "C", "C", "P1", "P2", "P3", "P4", "P5"),
    Intensities = c(150, 250, 225, 300, 175, 325, 600, 150, 350, 275),
    type = c(
      "Control", "Control", "Control", "Control", "Control",
      "Patients", "Patients", "Patients", "Patients", "Patients"
    ),
    group_size = c(5, 5, 5, 5, 5, 1, 1, 1, 1, 1)
  )

  test_hmdb_id <- "Test Metab 1"

  expect_silent(create_boxplot(test_ints_df_long, test_hmdb_id))

  expect_doppelganger(
    title = "create boxplot excel",
    fig = create_boxplot(test_ints_df_long, test_hmdb_id)
  )
})

testthat::test_that("save_plot_to_excel_workbook: Make and save a boxplot of intensities to an Excel workbook", {
  test_excel_workbook <- openxlsx::createWorkbook("Test")
  openxlsx::addWorksheet(test_excel_workbook, "Test_with_plots")

  test_intensities_df_long <- data.frame(
    Samples = c("C", "C", "C", "C", "C", "P1", "P2", "P3", "P4", "P5"),
    Intensities = c(150, 250, 225, 300, 175, 325, 600, 150, 350, 275),
    type = c(
      "Control", "Control", "Control", "Control", "Control",
      "Patients", "Patients", "Patients", "Patients", "Patients"
    ),
    group_size = c(5, 5, 5, 5, 5, 1, 1, 1, 1, 1)
  )

  test_hmdb_id <- "Test_metab_1"
  test_sheetname <- "Test_with_plots"
  test_file_path <- "./plot_test_"
  test_plot_width <- length(unique(test_intensities_df_long$Samples)) * 40
  test_col_width <- test_plot_width * 2
  test_start_row_index <- 1

  expect_silent(save_plot_to_excel_workbook(test_excel_workbook,
                                            test_sheetname,
                                            test_intensities_df_long,
                                            test_file_path,
                                            test_hmdb_id,
                                            test_plot_width,
                                            test_col_width,
                                            test_start_row_index))

  expect_identical(test_excel_workbook$media$image1.png, "./plot_test_Test_metab_1.png")

  # Remove test png file
  unlink(paste0(test_file_path, test_hmdb_id, ".png"))
  rm(test_excel_workbook)
})
