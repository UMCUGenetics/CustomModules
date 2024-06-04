## adapted from ident.hires.noise.HPC
# refactor: remove variables slope, incpt, ppm_iso_fixed
# combine with function get_element_info
# modified identify function to also look for adducts and their isotopes
identify_noisepeaks <- function(peakgroup_list, all_adducts, scanmode = "Negative", look4 = c("Cl", "Ac"),
                                noise_mz = NULL, resol = 140000, slope = 0, incpt = 0, ppm_fixed, ppm_iso_fixed) {
  #' Replace intensities that are zero with random value
  #'
  #' @param peakgroup_list: Peak group list (matrix)
  #' @param all_adducts: List of adducts to take into account (list of strings)
  #' @param scanmode: Scan mode, positive or negative (string)
  #' @param look4: List of adducts to look for (list of strings)
  #' @param noise_mz: All known noise peaks (matrix)
  #' @param resol: Value for resolution (integer)
  #' @param slope: Value for slope for mass correction (float)
  #' @param incpt: Value for intercept for mass correction (float)
  #' @param ppm_fixed: Value for distance between two values of mass (integer)
  #' @param ppm_iso_fixed: Value for distance between two values of mass for isotope peaks (integer)
  #'
  #' @return final_outlist: peak group list with filled-in intensities (matrix)

  options(stringsAsFactors = FALSE)
  metlin <- assi <- iso <- rep("", nrow(peakgroup_list))
  theormz <- nisos <- expint <- conf <- rep(0, nrow(peakgroup_list))

  # add adducts to identification list
  if (scanmode == "Positive") {
    adduct_scanmode <-  "+"
  } else {
    adduct_scanmode <- "-"
  }
  # make a copy of noise_mz
  noise_mz_orig <- noise_mz

  # loop over type of adduct
  for (adduct_index in 1:length(look4)) {
    noise_mz_adduct <- noise_mz_orig
    noise_mz_adduct[, "CompoundName"] <- as.character(noise_mz_orig[, "CompoundName"])

    if (look4[adduct_index] == "H2O") {
      add2label <- paste0("[M-", look4[adduct_index], "]", adduct_scanmode)
    } else {
      add2label <- paste0("[M+", look4[adduct_index], "]", adduct_scanmode)
    }

    noise_mz_adduct[, "CompoundName"] <- paste0(noise_mz_adduct[, "CompoundName"], add2label)
    adduct_info <- get_element_info(look4[adduct_index], all_adducts)
    if (scanmode == "Positive") {
      adduct_mass <- adduct_info$mass[1] + adduct_info$isotope$mass[1] - hydrogen_mass
    } else {
      adduct_mass <- adduct_info$mass[1] + adduct_info$isotope$mass[1] + hydrogen_mass
    }

    # loop over compounds in database
    for (compound_index in 1:nrow(noise_mz_adduct)) {
      # construct information for compound + adduct:
      if (scanmode == "Positive") {
        noise_mz_adduct[compound_index, "Mpos"] <-  as.numeric(noise_mz_adduct[compound_index, "Mpos"]) + adduct_mass
        noise_mz_adduct[compound_index, "MNeg"] <-  0
      } else {
        noise_mz_adduct[compound_index, "Mpos"] <-  0
        noise_mz_adduct[compound_index, "MNeg"] <-  as.numeric(noise_mz_adduct[compound_index, "MNeg"]) + adduct_mass
      }
    }
    noise_mz <- rbind(noise_mz, noise_mz_adduct)
  }

  if (scanmode == "Positive") {
    theor_mcol <- as.numeric(noise_mz[, "Mpos"])
  } else {
    theor_mcol <- as.numeric(noise_mz[, "MNeg"])
  }

  # get mz information from peakgroup_list
  mcol <- peakgroup_list[, "mzmed.pgrp"]
  # if column with average intensities is missing, calculate it:
  if (!("avg.int" %in% colnames(peakgroup_list))) {
    mzmaxcol <- which(colnames(peakgroup_list) == "mzmax.pgrp")
    endcol <- ncol(peakgroup_list)
    peakgroup_list[, "avg.int"] <- apply(peakgroup_list[, (mzmaxcol + 1):(endcol)], 1, mean)
  }

  # do indentification using own database:
  for (row_index in 1:nrow(noise_mz)) {
    theor_mz <- theor_mcol[row_index]

    # set tolerance for mz accuracy of main peak
    mtol <- theor_mz * ppm_fixed / 1000000
    # find main peak
    selp <- which(mcol > (theor_mz - mtol) & mcol < (theor_mz + mtol))
    # if there is more than one candidate peak for main, select best one based on mz_diff
    if (length(selp) > 1) {
      selp <- selp[abs(mcol[selp] - theor_mz) == min(abs(mcol[selp] - theor_mz))]
    }
    if (length(selp) == 1) {
      assi[selp] <- paste(assi[selp], as.character(noise_mz[row_index, "CompoundName"]), sep = ";")
      theormz[selp] <- theor_mz
    }
  }

  final_outlist <- cbind(peakgroup_list, assi, theormz, conf, nisos, iso, expint, metlin)
  return(final_outlist)
}
