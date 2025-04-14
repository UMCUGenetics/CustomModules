# unit tests for GenerateExcel
# functions: get_intensities_cols, calculate_zscores,
#            robust_scaler, remove_outliers_grubbs,
#            save_to_rdata_and_txt, set_row_height_col_width_wb
source("../../export/generate_excel_functions.R")

testthat::test_that("Get indices of columns and dataframe of intensities of a given label", {
  test_outlist <- data.frame(
    plots = c(NA, NA, NA, NA),
    C101.1 = c(100, 200, 300, 400),
    C102.1 = c(125, 225, 325, 425),
    P2.1 = c(150, 250, 350, 450),
    P3.1 = c(175, 275, 375, 475),
    HMDB_name = c("metab_1", "metab_2", "metab_3", "metab_4"),
    HMDB_name_all = c("metab_1;metab_11", "metab_2", "metab_3;metab_13", "metab_4"),
    HMDB_ID_all = c("HMDB001;HMDB011", "HMDB002", "HMDB003;HMDB013", "HMDB004"),
    sec_HMDB_ID = c("HMDB1;HMDB11", "", "HMDB3;HMDB13", "HMDB4"),
    HMDB_key = c("HMDB001", "HMDB002", "HMDB003", "HMDB004"),
    sec_HMDB_ID_rlvc = c("HMDB1 | HMDB11", "HMDB2", "HMDB3", "HMDB4"),
    name = c("metab_1 | metab_11", "metab_2", "metab_3", "metab_4"),
    relevance = c("Endogenous, relevant", "Endogenous, relevant | Exogenous", "Endogenous, relevant", "Endogenous, relevant"),
    descr = c("descr1", "descr2", "descr3", "descr4"),
    origin = c("Endogenous", "Endogenous | Exogenous", "Endogenous", "Endogenous"),
    fluids = c("Blood", "Blood", "Blood", "Blood"),
    tissue = c("Muscle", "", "Prostate", ""),
    disease = c("disease 1", "disease2", "", "disease3"),
    pathway = c("pathway1", "pathway2", "pathway3", "pathway4"),
    HMDB_code = c("HMDB001", "HMDB002", "HMDB003", "HMDB004")
  )
  rownames(test_outlist) <- c("HMDB001", "HMDB002", "HMDB003", "HMDB004")
  control_label <- "C"
  case_label <- "P"

  expect_equal(get_intensities_cols(test_outlist, control_label)$col_idx, c(2, 3))
  expect_equal(colnames(get_intensities_cols(test_outlist, control_label)$df_intensities), c("C101.1", "C102.1"))
  expect_equal(rownames(get_intensities_cols(test_outlist, control_label)$df_intensities), c("HMDB001", "HMDB002", "HMDB003", "HMDB004"))
  expect_equal(get_intensities_cols(test_outlist, control_label)$df_intensities$C101.1, c(100, 200, 300, 400))

  expect_equal(get_intensities_cols(test_outlist, case_label)$col_idx, c(4, 5))
  expect_equal(colnames(get_intensities_cols(test_outlist, case_label)$df_intensities), c("P2.1", "P3.1"))
})

testthat::test_that("Calculating Z-scores using different methods for excluding controls", {
  test_outlist <- data.frame(
    plots = c(NA, NA, NA, NA),
    C101.1 = c(1000, 1200, 1300, 1400),
    C102.1 = c(1100, 1700, 925, 1125),
    C103.1 = c(1300, 750, 1000, 1220),
    C104.1 = c(1650, 925, 1600, 1650),
    C105.1 = c(180000, 1950, 750, 15050),
    C106.1 = c(1050, 1100, 1200, 1300),
    C107.1 = c(1150, 1250, 825, 1025),
    C108.1 = c(1350, 850, 1175, 1420),
    C109.1 = c(1450, 1025, 1500, 1550),
    C110.1 = c(1200, 950, 1750, 19050),
    C111.1 = c(1050, 1125, 1300, 1450),
    C112.1 = c(1250, 975, 1450, 1950),
    P2.1 = c(3000, 12500, 5500, 6500),
    P3.1 = c(5750, 2750, 6750, 9750),
    HMDB_name = c("metab_1", "metab_2", "metab_3", "metab_4"),
    HMDB_name_all = c("metab_1;metab_11", "metab_2", "metab_3;metab_13", "metab_4"),
    HMDB_ID_all = c("HMDB001;HMDB011", "HMDB002", "HMDB003;HMDB013", "HMDB004"),
    sec_HMDB_ID = c("HMDB1;HMDB11", "", "HMDB3;HMDB13", "HMDB4"),
    HMDB_key = c("HMDB001", "HMDB002", "HMDB003", "HMDB004"),
    sec_HMDB_ID_rlvc = c(c("HMDB1 | HMDB11", "HMDB2", "HMDB3", "HMDB4")),
    name = c("metab_1 | metab_11", "metab_2", "metab_3", "metab_4"),
    relevance = c("Endogenous, relevant", "Endogenous, relevant | Exogenous", "Endogenous, relevant", "Endogenous, relevant"),
    descr = c("descr1", "descr2", "descr3", "descr4"),
    origin = c("Endogenous", "Endogenous | Exogenous", "Endogenous", "Endogenous"),
    fluids = c("Blood", "Blood", "Blood", "Blood"),
    tissue = c("Muscle", "", "Prostate", ""),
    disease = c("disease 1", "disease2", "", "disease3"),
    pathway = c("pathway1", "pathway2", "pathway3", "pathway4"),
    HMDB_code = c("HMDB001", "HMDB002", "HMDB003", "HMDB004")
  )
  rownames(test_outlist) <- c("HMDB001", "HMDB002", "HMDB003", "HMDB004")

  control_intensities <- data.frame(
    C101.1 = c(1000, 1200, 1300, 1400),
    C102.1 = c(1100, 1700, 925, 1125),
    C103.1 = c(1300, 750, 1000, 1220),
    C104.1 = c(1650, 925, 1600, 1650),
    C105.1 = c(180000, 1950, 750, 15050),
    C106.1 = c(1050, 1100, 1200, 1300),
    C107.1 = c(1150, 1250, 825, 1025),
    C108.1 = c(1350, 850, 1175, 1420),
    C109.1 = c(1450, 1025, 1500, 1550),
    C110.1 = c(1200, 950, 1750, 19050),
    C111.1 = c(1050, 1125, 1300, 1450),
    C112.1 = c(1250, 975, 1450, 1950)
  )
  rownames(control_intensities) <- c("HMDB001", "HMDB002", "HMDB003", "HMDB004")

  control_col_idx <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
  intensity_col_ids <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
  startcol <- ncol(test_outlist) + 3
  perc <- 5
  outlier_threshold <- 2

  expect_type(calculate_zscores(test_outlist, "_Zscore", control_intensities, NULL, intensity_col_ids, startcol), "list")
  expect_identical(colnames(calculate_zscores(test_outlist, "_Zscore", control_intensities, NULL, intensity_col_ids, startcol)),
                   c("plots", "C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "C106.1", "C107.1", "C108.1", "C109.1", "C110.1",
                     "C111.1", "C112.1", "P2.1", "P3.1", "HMDB_name", "HMDB_name_all", "HMDB_ID_all", "sec_HMDB_ID",
                     "HMDB_key", "sec_HMDB_ID_rlvc", "name", "relevance", "descr", "origin", "fluids", "tissue", "disease",
                     "pathway", "HMDB_code", "avg.ctrls", "sd.ctrls", "C101.1_Zscore", "C102.1_Zscore", "C103.1_Zscore",
                     "C104.1_Zscore", "C105.1_Zscore", "C106.1_Zscore", "C107.1_Zscore", "C108.1_Zscore", "C109.1_Zscore",
                     "C110.1_Zscore", "C111.1_Zscore", "C112.1_Zscore", "P2.1_Zscore", "P3.1_Zscore"))
  expect_equal(round(calculate_zscores(test_outlist, "_Zscore", control_intensities, NULL, intensity_col_ids, startcol)$avg.ctrls, 3),
               c(16129.167, 1150.0, 1231.250, 4015.833), tolerance = 0.001)
  expect_equal(calculate_zscores(test_outlist, "_Zscore", control_intensities, NULL, intensity_col_ids, startcol)$P2.1_Zscore,
               c(-0.2544103, 32.4586955, 13.6066674,  0.4037668), tolerance = 0.001)

  expect_type(calculate_zscores(test_outlist, "_RobustZscore", control_col_idx, perc, intensity_col_ids, startcol), "list")
  expect_identical(colnames(calculate_zscores(test_outlist, "_RobustZscore", control_col_idx, perc, intensity_col_ids, startcol)),
                   c("plots", "C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "C106.1", "C107.1", "C108.1", "C109.1", "C110.1",
                     "C111.1", "C112.1", "P2.1", "P3.1", "HMDB_name", "HMDB_name_all", "HMDB_ID_all", "sec_HMDB_ID",
                     "HMDB_key", "sec_HMDB_ID_rlvc", "name", "relevance", "descr", "origin", "fluids", "tissue", "disease",
                     "pathway", "HMDB_code", "avg.ctrls", "sd.ctrls", "C101.1_RobustZscore", "C102.1_RobustZscore",
                     "C103.1_RobustZscore", "C104.1_RobustZscore", "C105.1_RobustZscore", "C106.1_RobustZscore",
                     "C107.1_RobustZscore", "C108.1_RobustZscore", "C109.1_RobustZscore", "C110.1_RobustZscore",
                     "C111.1_RobustZscore", "C112.1_RobustZscore", "P2.1_RobustZscore", "P3.1_RobustZscore"))
  expect_equal(calculate_zscores(test_outlist, "_RobustZscore", control_col_idx, perc, intensity_col_ids, startcol)$avg.ctrls,
               c(1255.0, 1110.0, 1227.5, 2811.5), tolerance = 0.001)
  expect_equal(calculate_zscores(test_outlist, "_RobustZscore", control_col_idx, perc, intensity_col_ids, startcol)$P2.1_RobustZscore,
               c(9.1511750, 46.9804468, 16.8039663, 0.8565111), tolerance = 0.001)

  expect_type(calculate_zscores(test_outlist, "_OutlierRemovedZscore", control_col_idx, outlier_threshold, intensity_col_ids, startcol), "list")
  expect_identical(colnames(calculate_zscores(test_outlist, "_OutlierRemovedZscore", control_col_idx, outlier_threshold, intensity_col_ids, startcol)),
                   c("plots", "C101.1", "C102.1", "C103.1", "C104.1", "C105.1", "C106.1", "C107.1", "C108.1", "C109.1", "C110.1",
                     "C111.1", "C112.1", "P2.1", "P3.1", "HMDB_name", "HMDB_name_all", "HMDB_ID_all", "sec_HMDB_ID",
                     "HMDB_key", "sec_HMDB_ID_rlvc", "name", "relevance", "descr", "origin", "fluids", "tissue", "disease",
                     "pathway", "HMDB_code", "avg.ctrls", "sd.ctrls", "nr.ctrls", "C101.1_OutlierRemovedZscore",
                     "C102.1_OutlierRemovedZscore", "C103.1_OutlierRemovedZscore", "C104.1_OutlierRemovedZscore",
                     "C105.1_OutlierRemovedZscore", "C106.1_OutlierRemovedZscore", "C107.1_OutlierRemovedZscore",
                     "C108.1_OutlierRemovedZscore", "C109.1_OutlierRemovedZscore", "C110.1_OutlierRemovedZscore",
                     "C111.1_OutlierRemovedZscore", "C112.1_OutlierRemovedZscore", "P2.1_OutlierRemovedZscore",
                     "P3.1_OutlierRemovedZscore")
                    )
  expect_equal(calculate_zscores(test_outlist, "_OutlierRemovedZscore", control_col_idx, outlier_threshold, intensity_col_ids, startcol)$avg.ctrls,
               c(1231.818, 1077.273, 1231.250, 2649.091), tolerance = 0.001)
  expect_identical(calculate_zscores(test_outlist, "_OutlierRemovedZscore", control_col_idx, outlier_threshold, intensity_col_ids, startcol)$nr.ctrls,
                   c(11, 11, 12, 11))
  expect_equal(calculate_zscores(test_outlist, "_OutlierRemovedZscore", control_col_idx, outlier_threshold, intensity_col_ids, startcol)$P2.1_OutlierRemovedZscore,
               c(8.9955723, 44.9136860, 13.6066674, 0.9345077), tolerance = 0.001)
})

testthat::test_that("Use robust scaler", {
  control_intensities <- data.frame(
    C101.1 = c(1000),
    C102.1 = c(1100),
    C103.1 = c(1300),
    C104.1 = c(1650),
    C105.1 = c(180000),
    C106.1 = c(1050),
    C107.1 = c(1150),
    C108.1 = c(1350),
    C109.1 = c(1450),
    C110.1 = c(1200),
    C111.1 = c(1050),
    C112.1 = c(1250)
  )
  rownames(control_intensities) <- "HMDB001"
  control_col_idx <- c(1)
  perc <- 5

  expect_type(robust_scaler(control_intensities[1, 1:12], control_col_idx, perc), "double")
  expect_length(robust_scaler(control_intensities[1, 1:12], control_col_idx, perc), 10)
  expect_equal(robust_scaler(control_intensities[1, 1:12], control_col_idx, perc),
               c(1050, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1450, 1650), tolerance = 0.001)
})

testthat::test_that("Use Grubbs outlier removal", {
  control_intensities <- data.frame(
    C101.1 = c(1000),
    C102.1 = c(1100),
    C103.1 = c(1300),
    C104.1 = c(1650),
    C105.1 = c(180000),
    C106.1 = c(1050),
    C107.1 = c(1150),
    C108.1 = c(1350),
    C109.1 = c(1450),
    C110.1 = c(1200),
    C111.1 = c(1050),
    C112.1 = c(1250)
  )
  rownames(control_intensities) <- "HMDB001"
  outlier_threshold <- 2

  expect_type(remove_outliers_grubbs(control_intensities[1, ], outlier_threshold), "list")
  expect_length(remove_outliers_grubbs(control_intensities[1, ], outlier_threshold), 11)
  expect_equal(as.numeric(remove_outliers_grubbs(control_intensities[1, ], outlier_threshold)),
               c(1000, 1100, 1300, 1650, 1050, 1150, 1350, 1450, 1200, 1050, 1250), tolerance = 0.001)
  expect_identical(colnames(remove_outliers_grubbs(control_intensities[1, ], outlier_threshold)),
                   c("C101.1", "C102.1", "C103.1", "C104.1", "C106.1", "C107.1", "C108.1",
                     "C109.1", "C110.1", "C111.1", "C112.1"))
})

testthat::test_that("Save data to RData and txt file", {
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
  expect_true(file.exists("test_df.txt"))
  expect_true(file.exists("test_df.RData"))

  test_df_txt <- read.delim("test_df.txt")
  expect_identical(colnames(test_df_txt), c("C101.1", "C102.1", "P2.1", "P3.1", "HMDB_name", "sec_HMDB_ID"))
  expect_equal(test_df_txt$P2.1, c(150, 250, 350, 450))

  test_df_rdata <- get(load("test_df.RData"))
  expect_identical(colnames(test_df_rdata), c("C101.1", "C102.1", "P2.1", "P3.1", "HMDB_name", "sec_HMDB_ID"))
  expect_equal(test_df_rdata$P2.1, c(150, 250, 350, 450))

  file.remove("test_df.txt", "test_df.RData")
})

testthat::test_that("Check row height and column width in a workbook", {
  test_wb_plots <- openxlsx::createWorkbook("Test")
  openxlsx::addWorksheet(test_wb_plots, "Test_with_plots")

  sheetname_with_plots <- "Test_with_plots"
  num_rows_df <- 5
  num_cols_df <- 5
  plot_width <- 100

  correct_col_widths <- c("5", "20", "20", "20", "20")
  names(correct_col_widths) <- c(1, 2, 3, 4, 5)
  attr(correct_col_widths, "hidden") <- c("0", "0", "0", "0", "0")
  expect_identical(set_row_height_col_width_wb(test_wb_plots, sheetname_with_plots, num_rows_df, num_cols_df, plot_width,
                                               plots_present = TRUE)$colWidths[[1]], correct_col_widths)

  correct_row_heights <- c("140", "140", "140", "140", "140")
  names(correct_row_heights) <- c(2, 3, 4, 5, 6)
  expect_identical(set_row_height_col_width_wb(test_wb_plots, sheetname_with_plots, num_rows_df, num_cols_df, plot_width,
                                               plots_present = TRUE)$rowHeights[[1]], correct_row_heights)

  rm(test_wb_plots)

  test_wb_no_plots <- openxlsx::createWorkbook("Test")
  openxlsx::addWorksheet(test_wb_no_plots, "Test_no_plots")
  sheetname_no_plots <- "Test_no_plots"

  correct_col_widths <- c("20", "20", "20", "20", "20")
  names(correct_col_widths) <- c(1, 2, 3, 4, 5)
  attr(correct_col_widths, "hidden") <- c("0", "0", "0", "0", "0")
  expect_identical(set_row_height_col_width_wb(test_wb_no_plots, sheetname_no_plots, num_rows_df, num_cols_df, plot_width = NULL,
                                               plots_present = FALSE)$colWidths[[1]], correct_col_widths)

  correct_row_heights <- c("18", "18", "18", "18", "18")
  names(correct_row_heights) <- c(1, 2, 3, 4, 5)
  expect_identical(set_row_height_col_width_wb(test_wb_no_plots, sheetname_no_plots, num_rows_df, num_cols_df, plot_width = NULL,
                                               plots_present = FALSE)$rowHeights[[1]], correct_row_heights)

  rm(test_wb_no_plots)
})
