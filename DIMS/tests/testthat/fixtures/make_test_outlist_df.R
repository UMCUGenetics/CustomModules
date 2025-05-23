make_test_outlist_df <- function() {
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
  
  write.table(test_outlist, file = "test_outlist.txt", sep = "\t")
}

make_control_intensities_df <- function() {
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
  
  write.table(control_intensities, file = "tests/testthat/fixtures/test_control_intensities.txt", sep = "\t")
}

make_internal_standards_df <- function() {
  test_internal_standards <- data.frame(
    plots = c(NA, NA, NA, NA),
    C101.1 = c(100, 200, 300, 400),
    C102.1 = c(125, 225, 325, 425),
    P2.1 = c(150, 250, 350, 450),
    P3.1 = c(175, 275, 375, 475),
    HMDB_name = c("metab_1 (IS)", "metab_2 (IS)", "metab_3 (IS)", "metab_4 (IS)"),
    HMDB_name_all = c("metab_1 (IS)", "metab_2 (IS)", "metab_3 (IS)", "metab_4 (IS)"),
    HMDB_ID_all = c("HMDB100000", "HMDB200000", "HMDB300000", "HMDB400000"),
    sec_HMDB_ID = c("HMDB100000", "HMDB200000", "HMDB300000", "HMDB400000"),
    HMDB_key = c("HMDB100000", "HMDB200000", "HMDB300000", "HMDB400000"),
    sec_HMDB_ID_rlvc = c("HMDB100000", "HMDB200000", "HMDB300000", "HMDB400000"),
    name = c("metab_1 (IS)", "metab_2 (IS)", "metab_3 (IS)", "metab_4 (IS)"),
    relevance = c("Internal standard", "Internal standard", "Internal standard", "Internal standard"),
    descr = c("Internal standard", "Internal standard", "Internal standard", "Internal standard"),
    origin = c("Internalstandard", "Internalstandard", "Internalstandard", "Internalstandard"),
    fluids = c("Internalstandard", "Internalstandard", "Internalstandard", "Internalstandard"),
    tissue = c("Internalstandard", "Internalstandard", "Internalstandard", "Internalstandard"),
    disease = c("Internalstandard", "Internalstandard", "Internalstandard", "Internalstandard"),
    pathway = c("Internalstandard", "Internalstandard", "Internalstandard", "Internalstandard"),
    HMDB_code = c("HMDB100000", "HMDB200000", "HMDB300000", "HMDB400000"),
    avg.ctrls = c(137.5, 237.5, 337.5, 437.5),
    sd.ctrls = c(32.2749, 32.2749, 32.2749, 32.2749),
    nr.ctrls = c(2, 2, 2, 2),
    C101.1_Zscore = c(0.25, 1.65, 2.75, -0.35),
    C102.2_Zscore = c(1.89, -0.42, 0.22, 1.11),
    P2.1_Zscore = c(0.59, 1.36, -0.51, -0.32),
    P3.1_Zscore = c(1.03, 0.28, 0.78, 0.68)
  )
  rownames(test_internal_standards) <- c("HMDB100000", "HMDB200000", "HMDB300000", "HMDB400000")
  
  write.table(test_internal_standards, file = "tests/testthat/fixtures/test_internal_standards.txt", sep = "\t")
}

make_test_outlist_IS_df <- function() {
  test_outlist_IS <- data.frame(
    C101.1 = c(100, 200, 300, 400),
    C102.1 = c(125, 225, 325, 425),
    P2.1 = c(150, 250, 350, 450),
    P3.1 = c(175, 275, 375, 475),
    HMDB_name = c("metab_1 (IS)", "metab_85", "metab_3 (IS)", "metab_245"),
    HMDB_ID_all = c("HMDB100000", "HMDB68425", "HMDB300000", "HMDB84684"),
    sec_HMDB_ID = c("HMDB100000", "HMDB68425", "HMDB300000", "HMDB84684"),
    HMDB_name_all = c("metab_1 (IS)", "metab_85", "metab_3 (IS)", "metab_245")
  )
  rownames(test_outlist_IS) <- c("HMDB100000", "HMDB68425", "HMDB300000", "HMDB84684")
  
  write.table(test_outlist_IS, file = "tests/testthat/fixtures/test_outlist_IS.txt", sep = "\t")
}









