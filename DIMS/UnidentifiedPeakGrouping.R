#!/usr/bin/Rscript
## adapted from 8-peakGrouping.rest.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

unidentified_peaklist <- cmd_args[1]
resol <- as.numeric(cmd_args[2])
ppm <- as.numeric(cmd_args[3])
outdir <- "./"

options(digits = 16)

# function for grouping unidentified peaks
grouping_rest <- function(outdir, unidentified_peaklist, scanmode, ppm) {
  outlist_copy <- get(load(unidentified_peaklist))
  load(paste0(scanmode, "_repl_pattern.RData"))
  outpgrlist <- NULL

  # group on highest peaks
  range <- ppm * 1e-06

  # temporary: speed up this step by limiting the number of rows used in while loop
  # this script needs to be parallellized
  nrow_div <- nrow(outlist_copy) / 1.1

  # while (dim(outlist_copy)[1] > 0) {
  while (dim(outlist_copy)[1] > nrow_div) {
    sel <- which(as.numeric(outlist_copy[, "height.pkt"]) == max(as.numeric(outlist_copy[, "height.pkt"])))[1]

    # ppm range around max
    mzref <- as.numeric(outlist_copy[sel, "mzmed.pkt"])
    pkmin <- -(range * mzref - mzref)
    pkmax <- 2 * mzref - pkmin

    selp <- as.numeric(outlist_copy[, "mzmed.pkt"]) > pkmin & as.numeric(outlist_copy[, "mzmed.pkt"]) < pkmax
    tmplist <- outlist_copy[selp, , drop = FALSE]

    nrsamples <- length(unique(tmplist[, "samplenr"]))
    if (nrsamples > 0) {
      mzmed_pgrp <- mean(as.numeric(outlist_copy[selp, "mzmed.pkt"]))
      mzmin_pgrp <- -(range * mzmed_pgrp - mzmed_pgrp)
      mzmax_pgrp <- 2 * mzmed_pgrp - mzmin_pgrp
      # select peaks within mz range
      selp <- as.numeric(outlist_copy[, "mzmed.pkt"]) > mzmin_pgrp & as.numeric(outlist_copy[, "mzmed.pkt"]) < mzmax_pgrp
      tmplist <- outlist_copy[selp, , drop = FALSE]

      # remove used peaks
      tmp <- as.vector(which(tmplist[, "height.pkt"] == -1))
      if (length(tmp) > 0) tmplist <- tmplist[-tmp, , drop = FALSE]

      nrsamples <- length(unique(tmplist[, "samplenr"]))
      fq_worst_pgrp <- as.numeric(max(outlist_copy[selp, "fq"]))
      fq_best_pgrp  <- as.numeric(min(outlist_copy[selp, "fq"]))
      ints_allsamps <- rep(0, length(names(repl_pattern_filtered)))
      names(ints_allsamps) <- names(repl_pattern_filtered)

      # Check for each sample if multiple peaks exists, if so take the sum
      labels <- unique(tmplist[, "samplenr"])
      ints_allsamps[labels] <- as.vector(unlist(lapply(labels, function(x) {
        sum(as.numeric(tmplist[which(tmplist[, "samplenr"] == x), "height.pkt"]))
      })))

      # combine all information
      outpgrlist <- rbind(outpgrlist, c(mzmed_pgrp, fq_best_pgrp, fq_worst_pgrp, nrsamples, 
					mzmin_pgrp, mzmax_pgrp, ints_allsamps, NA, NA, NA, NA))
    }

    outlist_copy <- outlist_copy[-which(selp == TRUE), , drop = FALSE]
  }

  outpgrlist <- as.data.frame(outpgrlist)
  colnames(outpgrlist)[1:6] <- c("mzmed.pgrp", "fq.best", "fq.worst", "nrsamples", "mzmin.pgrp", "mzmax.pgrp")
  colnames(outpgrlist)[(length(repl_pattern_filtered) + 7):ncol(outpgrlist)] <- c("assi_HMDB", "iso_HMDB", 
										  "HMDB_code", "theormz_HMDB")

  return(outpgrlist)
}

# scanmodes <- c("positive", "negative")

# for (scanmode in scanmodes) {
#   # generate peak group lists of the unidentified peaks
#   unidentified_peaklist <- paste0("SpectrumPeaks_", scanmode, "_Unidentified.RData")
#   outpgrlist <- grouping_rest(outdir, unidentified_peaklist, scanmode, ppm = ppm)
#   write.table(outpgrlist, file = paste0("PeakGroupList_", scanmode, "_Unidentified.txt"))

#   # save output in RData format for further processing
#   save(outpgrlist, file=paste0("PeakGroupList_", scanmode, "_Unidentified.RData"))
# }

# determine appropriate scanmode based on unidentified_peaklist file
if (grepl("negative", basename(unidentified_peaklist))) {
  scanmode <- "negative"
} else if (grepl("positive", basename(unidentified_peaklist))) {
  scanmode <- "positive"
}

# generate peak group lists of the unidentified peaks
outpgrlist <- grouping_rest(outdir, unidentified_peaklist, scanmode, ppm = ppm)

# determine part number of unidentified_peaklist file
part_number <- gsub("\\D", "", basename(unidentified_peaklist))

# save output in txt format
write.table(outpgrlist, file = paste0("PeakGroupList_", scanmode, "_part_", part_number, "_Unidentified.txt"))

# save output in RData format for further processing
save(outpgrlist, file=paste0("PeakGroupList_", scanmode, "_part_", part_number, "_Unidentified.RData"))
