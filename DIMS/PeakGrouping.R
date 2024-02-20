#!/usr/bin/Rscript
# adapted from 6-peakGrouping.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

hmdb_part_file <- cmd_args[1]
ppm <- as.numeric(cmd_args[2])

options(digits = 16)

# load part of the HMDB
hmdb_add_iso <- get(load(hmdb_part_file))

# determine appropriate scanmode based on hmdb_part_file
if (grepl("negative", basename(hmdb_part_file))) {
  scanmode <- "negative"
} else if (grepl("positive", basename(hmdb_part_file))) {
  scanmode <- "positive"
}

# determine batch number of HMDB part file
batch_number <- strsplit(basename(hmdb_part_file), ".", fixed = TRUE)[[1]][2]

# load file with spectrum peaks
spec_peaks_file <- paste0("SpectrumPeaks_", scanmode, ".RData")
load(spec_peaks_file)
outlist_copy <- outlist_tot
rm(outlist_tot)

# load replication pattern
pattern_file <- paste0(scanmode, "_repl_pattern.RData")
load(pattern_file)

# determine appropriate column name in HMDB part
if (scanmode == "negative") {
  column_label <- "MNeg"
} else {
  column_label <- "Mpos"
}

# Initialize
peakgrouplist_identified <- NULL
list_of_peaks_used_in_peak_groups_identified <- NULL

# First find peak groups identified based on HMDB masses
while (dim(hmdb_add_iso)[1] > 0) {
  index <- 1

  # take one m/z value from the HMDB part and calculate mass tolerance
  reference_mass <- as.numeric(hmdb_add_iso[index, column_label])
  mass_tolerance <- (reference_mass * ppm) / 10^6

  # find the peaks in the dataset with corresponding m/z
  mzmed <- as.numeric(outlist_copy[, "mzmed.pkt"])
  selp <- which((mzmed > (reference_mass - mass_tolerance)) & (mzmed < (reference_mass + mass_tolerance)))
  tmplist <- outlist_copy[selp, , drop = FALSE]
  list_of_peaks_used_in_peak_groups_identified <- rbind(list_of_peaks_used_in_peak_groups_identified, tmplist)

  nrsamples <- length(selp)
  if (nrsamples > 0) {
    mzmed_pgrp <- mean(as.numeric(outlist_copy[selp, "mzmed.pkt"]))
    mzmin_pgrp <- reference_mass - mass_tolerance
    mzmax_pgrp <- reference_mass + mass_tolerance

    # determine fit quality fq
    fq_worst_pgrp <- as.numeric(max(outlist_copy[selp, "fq"]))
    fq_best_pgrp <- as.numeric(min(outlist_copy[selp, "fq"]))

    # set up object for intensities for all samples
    ints_allsamps <- rep(0, length(names(repl_pattern_filtered)))
    names(ints_allsamps) <- names(repl_pattern_filtered)

    # Check for each sample if multiple peaks exist, if so take the sum of the intensities
    labels <- unique(tmplist[, "samplenr"])
    ints_allsamps[labels] <- as.vector(unlist(lapply(labels, function(x) {
      sum(as.numeric(tmplist[which(tmplist[, "samplenr"] == x), "height.pkt"]))
    })))

    # Initialize
    assi_hmdb <- iso_hmdb <- hmdb_code <- NA
    tmplist_mass_iso <- tmplist_mass_adduct <- NULL

    # Identification: find all entries in HMDB part with mass within ppm range
    mass_all <- as.numeric(hmdb_add_iso[, column_label])
    index <- which((mass_all > (reference_mass - mass_tolerance)) & (mass_all < (reference_mass + mass_tolerance)))
    tmplist_mass <- hmdb_add_iso[index, , drop = FALSE]

    if (dim(tmplist_mass)[1] > 0) {
      # find isotope entries
      index_iso <- grep(" iso ", tmplist_mass[, "CompoundName"], fixed = TRUE)
      if (length(index_iso) > 0) {
        tmplist_mass_iso <- tmplist_mass[index_iso, , drop = FALSE]
        tmplist_mass <- tmplist_mass[-index_iso, , drop = FALSE]
      }

      if (dim(tmplist_mass)[1] > 0) {
        # find adduct entries
        index_adduct <- grep(" [M", tmplist_mass[, "CompoundName"], fixed = TRUE)
        if (length(index_adduct) > 0) {
          tmplist_mass_adduct <- tmplist_mass[index_adduct, , drop = FALSE]
          tmplist_mass <- tmplist_mass[-index_adduct, , drop = FALSE]
        }
      }

      # Compose a list compounds, adducts or isotopes with corresponding m/z
      if (dim(tmplist_mass)[1] > 0) {
        # metabolites
        assi_hmdb <- as.character(paste(as.character(tmplist_mass[, "CompoundName"]), collapse = ";"))
        hmdb_code <- as.character(paste(as.character(rownames(tmplist_mass)), collapse = ";"))
        theormz_hmdb <- as.numeric(tmplist_mass[1, column_label])

        # adducts of metabolites
        if (!is.null(tmplist_mass_adduct)) {
          if (dim(tmplist_mass_adduct)[1] > 0) {
            if (is.na(assi_hmdb)) {
              assi_hmdb <- as.character(paste(as.character(tmplist_mass_adduct[, "CompoundName"]), collapse = ";"))
              hmdb_code <- as.character(paste(as.character(rownames(tmplist_mass_adduct)), collapse = ";"))
            } else {
              assi_hmdb <- paste(assi_hmdb,
                                 as.character(paste(as.character(tmplist_mass_adduct[, "CompoundName"]), collapse = ";")), sep = ";")
              hmdb_code <- paste(hmdb_code,
                                 as.character(paste(as.character(rownames(tmplist_mass_adduct)), collapse = ";")), sep = ";")
            }
          }
        }

        # isotopes of metabolites
        if (!is.null(tmplist_mass_iso)) {
          if (dim(tmplist_mass_iso)[1] > 0) {
            iso_hmdb <- as.character(paste(as.character(tmplist_mass_iso[, "CompoundName"]), collapse = ";"))
          }
        }

        # if no metabolites have the correct m/z, look for adducts and isotopes only
      } else if (!is.null(tmplist_mass_adduct)) {
        theormz_hmdb <- as.numeric(tmplist_mass_adduct[1, column_label])

        # adducts of metabolites
        if (!is.null(tmplist_mass_adduct)) {
          if (dim(tmplist_mass_adduct)[1] > 0) {
            if (is.na(assi_hmdb)) {
              assi_hmdb <- as.character(paste(as.character(tmplist_mass_adduct[, "CompoundName"]), collapse = ";"))
              hmdb_code <- as.character(paste(as.character(rownames(tmplist_mass_adduct)), collapse = ";"))
            } else {
              assi_hmdb <- paste(assi_hmdb,
                                 as.character(paste(as.character(tmplist_mass_adduct[, "CompoundName"]), collapse = ";")), sep = ";")
              hmdb_code <- paste(hmdb_code,
                                 as.character(paste(as.character(rownames(tmplist_mass_adduct)), collapse = ";")), sep = ";")
            }
          }
        }

        # isotopes of metabolites
        if (!is.null(tmplist_mass_iso)) {
          if (dim(tmplist_mass_iso)[1] > 0) {
            iso_hmdb <- as.character(paste(as.character(tmplist_mass_iso[, "CompoundName"]), collapse = ";"))
          }
        }

        # if no metabolites or adducts can be found, only look for isotopes
      } else if (!is.null(tmplist_mass_iso)) {
        if (dim(tmplist_mass_iso)[1] > 0) {
          theormz_hmdb <- as.numeric(tmplist_mass_iso[1, column_label])
          iso_hmdb <- as.character(paste(as.character(tmplist_mass_iso[, "CompoundName"]), collapse = ";"))
        }
      }
    }

    # combine all information
    peakgrouplist_identified <- rbind(peakgrouplist_identified, cbind(
      data.frame(mzmed_pgrp, "fq.best" = fq_best_pgrp, "fq.worst" = fq_worst_pgrp, nrsamples, mzmin_pgrp, mzmax_pgrp),
      t(as.matrix(ints_allsamps)),
      data.frame(assi_hmdb, iso_hmdb, hmdb_code, theormz_hmdb)
    ))
  }

  # remove index metabolite from HMDB part and continue while loop
  hmdb_add_iso <- hmdb_add_iso[-index, ]
}


# save peak list corresponding to masses in HMDB part
save(list_of_peaks_used_in_peak_groups_identified, file = paste0(batch_number, "_", scanmode, "_peaks_used.RData"))
# save peak group list, identified part
save(peakgrouplist_identified, file = paste0(batch_number, "_", scanmode, "_identified.RData"))
