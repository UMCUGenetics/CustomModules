## adapted from replaceZeros.R
# this function does two things: replace zeros with random value and identify noise peaks
# refactor: split into two functions
# remove parameters outdir and thresh
# make hard-coded path to file with noise peaks into variable
# remove variables outdir, thresh
replace_zeros <- function(peakgroup_list, repl_pattern, scanmode, resol, outdir, thresh, ppm) {
  #' Replace intensities that are zero with random value
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #' @param repl_pattern: Replication pattern (list of strings)
  #' @param scanmode: Scan mode, positive or negative (string)
  #' @param resol: Value for resolution (integer)
  #' @param outdir: Path for output directory (string)
  #' @param thresh: Value for threshold (integer)
  #' @param ppm: Value for distance between two values of mass (integer)
  #'
  #' @return final_outlist: peak group list with filled-in intensities (matrix)

  # replace zeros
  if (!is.null(peakgroup_list)) {
    for (sample_index in 1:length(names(repl_pattern))) {
      sample_peaks <- peakgroup_list[, names(repl_pattern)[sample_index]]
      zero_intensity <- which(sample_peaks <= 0)
      if (!length(zero_intensity)) {
        next
      }
      for (zero_index in 1:length(zero_intensity)) {
        area <- fit_optim(peakgroup_list[zero_intensity[zero_index], "mzmed.pgrp"], thresh,
                          resol, FALSE, scanmode, int_factor = 1 * 10^5, 1, 1)$area
        peakgroup_list[zero_intensity[zero_index], names(repl_pattern)[sample_index]] <- rnorm(n = 1, mean = area,
                                                                                               sd = 0.25 * area)
      }
    }

    # Add column with average intensity
    peakgroup_list <- cbind(peakgroup_list, "avg.int" = apply(peakgroup_list[, 7:(ncol(peakgroup_list) - 4)], 1, mean))

    if (scanmode == "negative") {
      label <- "MNeg"
      label2 <- "Negative"
      # look for adducts in negative mode
      look4_adducts <- c("Cl", "Cl37", "For", "NaCl", "KCl", "H2PO4", "HSO4", "Na-H", "K-H", "H2O", "I")
    } else {
      label <- "Mpos"
      label2 <- "Positive"
      # look for adducts in positive mode
      look4_adducts <- c("Na", "K", "NaCl", "NH4", "2Na-H", "CH3OH", "KCl", "NaK-H")
    }

    # Identify noise peaks
    noise_mz <- read.table(file = "/hpc/dbg_mz/tools/db/TheoreticalMZ_NegPos_incNaCl.txt",
                           sep = "\t", header = TRUE, quote = "")
    noise_mz <- noise_mz[(noise_mz[, label] != 0), 1:4]
    outlist_withnoise <- identify_noisepeaks(peakgroup_list, all_adducts, scanmode = label2,
                                             noise_mz, look4 = look4_adducts, resol = resol,
                                             slope = 0, incpt = 0, ppm_fixed = ppm, ppm_iso_fixed = ppm)
    noise_info <- outlist_withnoise[, c("assi", "theormz")]
    colnames(noise_info) <- c("assi_noise",  "theormz_noise")

    final_outlist <- cbind(peakgroup_list, noise_info)
    return(final_outlist)
  }
}
