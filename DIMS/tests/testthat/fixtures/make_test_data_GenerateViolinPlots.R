### Functions used to create mock dataframes used for unit testing of GenerateViolinPlots ###

make_outlist_df <- function() {
  test_outlist_df <- data.frame(
    plots = NA,
    C101.1 = c(1000, 1200, 1300, 1400, 1500, 1600),
    C102.1 = c(1100, 1700, 925, 1125, 1200, 1050),
    C103.1 = c(1300, 750, 1000, 1220, 1100, 1200),
    C104.1 = c(1650, 925, 1600, 1650, 1025, 1150),
    C105.1 = c(180000, 1950, 750, 15050, 1100, 1300),
    P2025M1 = c(1000, 1200, 1300, 1400, 1100, 975),
    P2025M2 = c(1100, 1700, 925, 1125, 1050, 1175),
    P2025M3 = c(1300, 750, 1000, 1220, 975, 1100),
    P2025M4 = c(1650, 925, 1600, 1650, 1700, 1750),
    P2025M5 = c(180000, 1950, 750, 15050, 10000, 1500),
    HMDB_name = c("metab1", "metab2", "metab3", "metab4", "metab5", "metab6"),
    HMDB_name_all = NA,
    HMDB_ID_all = NA,
    sec_HMDB_ID = NA,
    HMDB_key = NA,
    sec_HMDB_ID_rlvnc = NA,
    name = NA,
    relevance = NA,
    descr = NA,
    origin = NA,
    fluids = NA,
    tissue = NA,
    disease = NA,
    pathway = NA,
    HMDB_code = c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012"),
    avg_ctrls = c(37010, 1305, 1115, 4089, 1185, 1260),
    sd_ctrls = c(79934.23, 508.80, 336.15, 6130.65, 186.75, 210.36),
    nr_ctrls = c(25, 26, 27, 28, 29, 30),
    C101.1_Zscore = c(0.45, 1.67, -1.86, 0.58, 2.47, -0.56),
    C102.1_Zscore = c(2.89, 0.79, -1.88, 5.46, -0.68, 1.65),
    C103.1_Zscore = c(0.54, -0.85, 1.58, 3.84, 0.84, -1.11),
    C104.1_Zscore = c(0.53, 1.84, 0.35, -0.54, 1.48, 0.43),
    C105.1_Zscore = c(3.46, -1.31, 0.14, -0.15, 1.48, 0.36),
    P2025M1_Zscore = c(0.31, 1.84, 2.34, 0.84, -0.46, 0.14),
    P2025M2_Zscore = c(2.45, 0.48, 1.45, -0.15, -1.51, 3.56),
    P2025M3_Zscore = c(2.14, 0.15, -1.44, -0.78, 1.68, 0.51),
    P2025M4_Zscore = c(12.18, 2.48, -0.18, 0.84, 1.48, -2.45),
    P2025M5_Zscore = c(3.22, 0.48, -3.18, 0.47, 1.18, 2.14)
  )
  rownames(test_outlist_df) <- c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012")
  
  write.table(test_outlist_df, file = "tests/testthat/fixtures/GenerateViolinPlots/test_outlist_df.txt", sep = "\t")
}

make_intensities_zscore_df <- function() {
  test_intensities_zscore_df <- data.frame(
    HMDB_code = c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012"),
    HMDB_name = c("metab1", "metab2", "metab3", "metab4", "metab5", "metab6"),
    C101.1 = c(1000, 1200, 1300, 1400, 1500, 1600),
    C102.1 = c(1100, 1700, 925, 1125, 1200, 1050),
    C103.1 = c(1300, 750, 1000, 1220, 1100, 1200),
    C104.1 = c(1650, 925, 1600, 1650, 1025, 1150),
    C105.1 = c(180000, 1950, 750, 15050, 1100, 1300),
    P2025M1 = c(1000, 1200, 1300, 1400, 1100, 975),
    P2025M2 = c(1100, 1700, 925, 1125, 1050, 1175),
    P2025M3 = c(1300, 750, 1000, 1220, 975, 1100),
    P2025M4 = c(1650, 925, 1600, 1650, 1700, 1750),
    P2025M5 = c(180000, 1950, 750, 15050, 10000, 1500),
    mean_control = c(37010, 1305, 1115, 4089, 1185, 1260),
    sd_control = c(79934.23, 508.80, 336.15, 6130.65, 186.75, 210.36),
    C101.1_Zscore = c(0.45, 1.67, -1.86, 0.58, 2.47, -0.56),
    C102.1_Zscore = c(2.89, 0.79, -1.88, 5.46, -0.68, 1.65),
    C103.1_Zscore = c(0.54, -0.85, 1.58, 3.84, 0.84, -1.11),
    C104.1_Zscore = c(0.53, 1.84, 0.35, -0.54, 1.48, 0.43),
    C105.1_Zscore = c(3.46, -1.31, 0.14, -0.15, 1.48, 0.36),
    P2025M1_Zscore = c(0.31, 1.84, 2.34, 0.84, -0.46, 0.14),
    P2025M2_Zscore = c(2.45, 0.48, 1.45, -0.15, -1.51, 3.56),
    P2025M3_Zscore = c(2.14, 0.15, -1.44, -0.78, 1.68, 0.51),
    P2025M4_Zscore = c(12.18, 2.48, -0.18, 0.84, 1.48, -2.45),
    P2025M5_Zscore = c(3.22, 0.48, -3.18, 0.47, 1.18, 2.14)
  )
  rownames(test_intensities_zscore_df) <- c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012")

  write.table(test_intensities_zscore_df, file = "tests/testthat/fixtures/test_intensities_zscore_df.txt", sep = "\t")
}

make_files_metabolite_groups <- function() {
  # Create new directories
  dir.create("tests/testthat/fixtures/test_metabolite_groups")

  test_acyl_carnitines <- data.frame(
    HMDB_code = c("HMDB001", "HMDB003", "HMDBAA1"),
    HMDB_name = c("metab1", "metab3", "ratio1"),
    Helix = c("ja", "nee", "ja"),
    Helix_naam = c("Metab_1", "Metab_3", "Ratio_1"),
    high_zscore = c(2, 2, 2),
    low_zscore = c(-1.5, -1.5, -1.5)
  )

  write.table(test_acyl_carnitines, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_metabolite_groups/test_acyl_carnitines.txt")

  test_crea_gua <- data.frame(
    HMDB_code = c("HMDB004", "HMDB011"),
    HMDB_name = c("metab4", "metab11"),
    Helix = c("ja", "ja"),
    Helix_naam = c("Metab_4", "Metab_11"),
    high_zscore = c(2, 2),
    low_zscore = c(-1.5, -1.5)
  )

  write.table(test_crea_gua, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_metabolite_groups/test_crea_gua.txt")
}

make_test_metab_interest_sort <- function() {
  test_acyl_carnitines_patients <- data.frame(
    HMDB_name = c(
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       "
    ),
    Sample = c("P2025M1", "P2025M1", "P2025M2", "P2025M2", "P2025M3",
               "P2025M3", "P2025M4", "P2025M4", "P2025M5", "P2025M5"),
    Z_score = c(0.31, 2.34, 2.45, 1.45, 2.14, -1.44, 12.18, -0.18, 3.22, -3.18)
  )

  test_acyl_carnitines_controls <- data.frame(
    HMDB_name = c(
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       "
    ),
    Sample = c("C101.1", "C101.1", "C102.1", "C102.1", "C103.1", "C103.1", "C104.1", "C104.1", "C105.1", "C105.1"),
    Z_score = c(0.45, -1.86, 2.89, -1.88, 0.54, 1.58, 0.53, 0.35, 3.46, 0.14)
  )

  write.table(test_acyl_carnitines_patients, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_acyl_carnitines_patients.txt")
  write.table(test_acyl_carnitines_controls, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_acyl_carnitines_controls.txt")

  test_crea_gua_patients <- data.frame(
    HMDB_name = c(
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      "
    ),
    Sample = c("P2025M1", "P2025M1", "P2025M2", "P2025M2", "P2025M3",
               "P2025M3", "P2025M4", "P2025M4", "P2025M5", "P2025M5"),
    Z_score = c(0.84, -0.46, -0.15, -1.51, -0.78, 1.68, 0.84, 1.48, 0.47, 1.18)
  )

  test_crea_gua_controls <- data.frame(
    HMDB_name = c(
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      "
    ),
    Sample = c("C101.1", "C101.1", "C102.1", "C102.1", "C103.1", "C103.1", "C104.1", "C104.1", "C105.1", "C105.1"),
    Z_score = c(0.58, 2.47, 5.46, -0.68, 3.84, 0.84, -0.54, 1.48, -0.15, 1.48)
  )

  write.table(test_crea_gua_patients, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_crea_gua_patients.txt")
  write.table(test_crea_gua_controls, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_crea_gua_controls.txt")
}

make_test_df_metabs_helix <- function() {
  test_df_metabs_helix <- data.frame(
    HMDB_name = c("metab1", "metab1", "metab1", "metab1", "metab1", "metab4", "metab11", "metab4",
                  "metab11", "metab4", "metab11", "metab4", "metab11", "metab4", "metab11"),
    Sample = c("P2025M1", "P2025M2", "P2025M3", "P2025M4", "P2025M5", "P2025M1", "P2025M1", "P2025M2", "P2025M2", "P2025M3",
               "P2025M3", "P2025M4", "P2025M4", "P2025M5", "P2025M5"),
    Z_score = c(0.31, 2.45, 2.14, 12.18, 3.22, 0.84, -0.46, -0.15, -1.51, -0.78, 1.68, 0.84, 1.48, 0.47, 1.18),
    Helix_naam = c("Metab_1", "Metab_1", "Metab_1", "Metab_1", "Metab_1", "Metab_4", "Metab_11",
                   "Metab_4", "Metab_11", "Metab_4", "Metab_11", "Metab_4", "Metab_11", "Metab_4", "Metab_11"),
    high_zscore = rep(2, 15),
    low_zscore = rep(-1.5, 15)
  )

  write.table(test_df_metabs_helix, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_df_metabs_helix.txt")
}

make_test_zscore_patient_df <- function() {
  test_zscore_patient_df <- data.frame(
    HMDB_code = c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB005", "HMDB006", "HMDB007", "HMDB008", "HMDB009", "HMDB010",
                  "HMDB011", "HMDB012", "HMDB013", "HMDB014", "HMDB015", "HMDB016", "HMDB017", "HMDB018", "HMDB019", "HMDB020",
                  "HMDB021", "HMDB022", "HMDB023", "HMDB024", "HMDB025", "HMDB026", "HMDB027", "HMDB028", "HMDB029",
                  "HMDB030"),
    HMDB_name = c("metab1", "metab2", "metab3", "metab4", "metab5", "metab6", "metab7", "metab8", "metab9", "metab10",
                  "metab11", "metab12", "metab13", "metab14", "metab15", "metab16", "metab17", "metab18", "metab19", "metab20",
                  "metab21", "metab22", "metab23", "metab24", "metab25", "metab26", "metab27", "metab28", "metab29",
                  "metab30"),
    P2025M1 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
    P2025M2 = c(-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15,
                -16, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28, -29, -30),
    P2025M3 = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                -16, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28, -29, -30),
    P2025M4 = c(-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15,
                16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
  )

  write.table(test_zscore_patient_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_zscore_patient_df.txt")
}

make_test_metab_perpage <- function() {
  test_acyl_carnitines_df <- data.frame(
    HMDB_name = c(
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       ",
      "metab1                                       ",
      "metab3                                       "
    ),
    Sample = c("P2025M1", "P2025M1", "P2025M2", "P2025M2", "P2025M3",
               "P2025M3", "P2025M4", "P2025M4", "P2025M5", "P2025M5",
               "C101.1", "C101.1", "C102.1", "C102.1", "C103.1", "C103.1",
               "C104.1", "C104.1", "C105.1", "C105.1"),
    Z_score = c(0.31, 2.34, 2.45, 1.45, 2.14, -1.44, 12.18, -0.18, 3.22, -3.18,
                0.45, -1.86, 2.89, -1.88, 0.54, 1.58, 0.53, 0.35, 3.46, 0.14)
  )

  write.table(test_acyl_carnitines_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_acyl_carnitines_df.txt")

  test_crea_gua_df <- data.frame(
    HMDB_name = c(
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      ",
      "metab4                                       ",
      "metab11                                      "
    ),
    Sample = c("P2025M1", "P2025M1", "P2025M2", "P2025M2", "P2025M3",
               "P2025M3", "P2025M4", "P2025M4", "P2025M5", "P2025M5",
               "C101.1", "C101.1", "C102.1", "C102.1", "C103.1", "C103.1",
               "C104.1", "C104.1", "C105.1", "C105.1"),
    Z_score = c(0.84, -0.46, -0.15, -1.51, -0.78, 1.68, 0.84, 1.48, 0.47, 1.18,
                0.58, 2.47, 5.46, -0.68, 3.84, 0.84, -0.54, 1.48, -0.15, 1.48)
  )

  write.table(test_crea_gua_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_crea_gua_df.txt")
}

make_test_expected_biomark_df <- function() {
  test_expected_biomarkers_df <- data.frame(
    HMDB_code = c("HMDB002", "HMDB002", "HMDB005", "HMDB005", "HMDB005", "HMDB009", "HMDB012", "HMDB012", "HMDB020", "HMDB020",
                  "HMDB025", "HMDB025", "HMDB025", "HMDB028", "HMDB028"),
    HMDB_name = c("metab2", "metab2", "metab5", "metab5", "metab5", "metab9", "metab12", "metab12", "metab20", "metab20",
                  "metab25", "metab25", "metab25", "metab28", "metab28"),
    Disease = c("Disease A", "Disease B", "Disease B", "Disease B", "Disease C", "Disease D", "Disease A", "Disease E",
                "Disease F", "Disease C", "Disease F", "Disease G", "Disease D", "Disease G", "Disease E"),
    M.z = c("1.2", "1.2", "2.5", "2.5", "2.5", "3.0", "4.5", "4.5", "2.5", "2.5", "5.1", "5.1", "5.1", "5.8", "5.8"),
    Change = c("Increase", "Decrease", "Increase", "Increase", "Decrease", "Decrease", "Increase", "Increase", "Decrease",
               "Increase", "Increase", "Decrease", "Increase", "Decrease", "Increase"),
    Total_Weight = c(10.0, -1.5, 2.0, 5.0, -2.5, -3.0, 14.5, 4.0,
                     -7.5, 6.0, 20.0, -5.0, 3.0, -1.5, 4.0),
    Absolute_Weight = c(10.0, 1.5, 2.0, 5.0, 2.5, 3.0, 14.5, 4.0,
                        7.5, 6.0, 20.0, 5.0, 3.0, 1.5, 4.0),
    Dispensability = c("Dispensable", "Dispensable", "Indispensable", "Dispensable", "Dispensable", "Indispensable",
                       "Dispensable", "Dispensable", "Indispensable", "Indispensable", "Dispensable", "Dispensable",
                       "Dispensable", "Dispensable", "Indispensable")
  )

  write.table(test_expected_biomarkers_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_expected_biomarkers_df.txt")
}

make_test_probability_score_df <- function() {
  test_probability_score_df <- data.frame(
    Disease = c("Disease A", "Disease B", "Disease C", "Disease D", "Disease E", "Disease F", "Disease G"),
    P2025M1 = c(10.9, 0.953, 12.1, 0, 44.3, 0, -38.7),
    P2025M2 = c(-10.9, 0, 0, -12.5, 0, -77.4, 38.7),
    P2025M3 = c(49.9, 2.29, 0, 0, 0, -77.4, 38.7),
    P2025M4 = c(-49.9, 0, 12.1, 18.2, 28.1, 0, -38.7)
  )

  write.table(test_probability_score_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/test_probability_score_df.txt")
}

make_zscore_dfs <- function() {
  test_zscore_patients_df <- data.frame(
    HMDB_code = c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3"),
    HMDB_name = c("metab1", "metab2", "metab3", "metab4", "metab5", "metab6", "Test_ratio1", "Test_ratio2", "Test_ratio3"),
    P2025M1 = c(0.31, 1.84, 2.34, 0.84, -0.46, 0.14, -0.58, 0.48, -0.45),
    P2025M2 = c(2.45, 0.48, 1.45, -0.15, -1.51, 3.56, -0.71, 0.39, -0.54),
    P2025M3 = c(2.14, 0.15, -1.44, -0.78, 1.68, 0.51, -0.22, 0.38, -0.48),
    P2025M4 = c(12.18, 2.48, -0.18, 0.84, 1.48, -2.45, -0.21, 0.51, -0.29),
    P2025M5 = c(3.22, 0.48, -3.18, 0.47, 1.18, 2.14, 1.74, -1.78, 1.78)
  )
  
  test_zscore_controls_df <- data.frame(
    HMDB_code = c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3"),
    HMDB_name = c("metab1", "metab2", "metab3", "metab4", "metab5", "metab6", "Test_ratio1", "Test_ratio2", "Test_ratio3"),
    C101.1 = c(0.45, 1.67, -1.86, 0.58, 2.47, -0.56, -0.58, 0.48, -0.45),
    C102.1 = c(2.89, 0.79, -1.88, 5.46, -0.68, 1.65, -0.71, 0.39, -0.54),
    C103.1 = c(0.54, -0.85, 1.58, 3.84, 0.84, -1.11, -0.22, 0.38, -0.48),
    C104.1 = c(0.53, 1.84, 0.35, -0.54, 1.48, 0.43, -0.21, 0.51, -0.29),
    C105.1 = c(3.46, -1.31, 0.14, -0.15, 1.48, 0.36, 1.74, -1.78, 1.78)
  )
  
  write.table(test_zscore_patients_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/GenerateViolinPlots/test_zscore_patients_df.txt")
  write.table(test_zscore_controls_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/GenerateViolinPlots/test_zscore_controls_df.txt")
}

make_test_metabolite_class_dfs <- function() {
  test_metabolite_class_patients_df <- data.frame(
    HMDB_name = c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3",
                  "HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3",
                  "HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3"),
    Sample = c("P2025M1", "P2025M1", "P2025M1", "P2025M1", "P2025M1", "P2025M1", "P2025M1", "P2025M1", "P2025M1",
               "P2025M2", "P2025M2", "P2025M2", "P2025M2", "P2025M2", "P2025M2", "P2025M2", "P2025M2", "P2025M2",
               "P2025M3", "P2025M3", "P2025M3", "P2025M3", "P2025M3", "P2025M3", "P2025M3", "P2025M3", "P2025M3"),
    Z_score = c(0.31, 1.84, 2.34, 0.84, -0.46, 0.14, -0.58, 0.48, -0.45,
                2.45, 0.48, 1.45, -0.15, -1.51, 3.56, -0.71, 0.39, -0.54,
                2.14, 0.15, -1.44, -0.78, 1.68, 0.51, -0.22, 0.38, -0.48)
  )
  
  test_metabolite_class_controls_df <- data.frame(
    HMDB_name = c("HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3",
                  "HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3",
                  "HMDB001", "HMDB002", "HMDB003", "HMDB004", "HMDB011", "HMDB012", "HMDB000TT1", "HMDB000TT2", "HMDB000TT3"),
    Sample = c("C101.1", "C101.1", "C101.1", "C101.1", "C101.1", "C101.1", "C101.1", "C101.1", "C101.1",
               "C102.1", "C102.1", "C102.1", "C102.1", "C102.1", "C102.1", "C102.1", "C102.1", "C102.1",
               "C103.1", "C103.1", "C103.1", "C103.1", "C103.1", "C103.1", "C103.1", "C103.1", "C103.1"),
    Z_score = c(0.45, 1.67, -1.86, 0.58, 2.47, -0.56, -0.58, 0.48, -0.45,
                2.89, 0.79, -1.88, 5.46, -0.68, 1.65, -0.71, 0.39, -0.54,
                0.54, -0.85, 1.58, 3.84, 0.84, -1.11, -0.22, 0.38, -0.48)
  )
  
  write.table(test_metabolite_class_patients_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/GenerateViolinPlots/test_metabolite_class_patients_df.txt")
  write.table(test_metabolite_class_controls_df, sep = "\t", quote = FALSE, row.names = FALSE,
              file = "tests/testthat/fixtures/GenerateViolinPlots/test_metabolite_class_controls_df.txt")
  
}
