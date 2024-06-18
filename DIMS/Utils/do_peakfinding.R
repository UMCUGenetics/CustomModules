## adapted from findpeaks.Gauss.HPC.R
# NB: this function will be taken up into PeakFinding.R
# variables with fixed values will be removed from function parameters
# int_factor, scale, outdir, plot, thresh, width, height
do_peakfinding <- function(sample_avgtechrepl, int_factor, scale, resol, outdir, scanmode, plot, thresh, width, height) {
  #' start peak finding
  #'
  #' @param sample_avgtechrepl: Dataframe with binned intensities averaged over technical replicates for a sample
  #' @param int_factor: Value used to calculate area under Gaussian curve (integer)
  #' @param scale: Initial value used to estimate scaling parameter (integer)
  #' @param resol: Value for resolution (integer)
  #' @param outdir: Path for output directory (string)
  #' @param scanmode: Scan mode, positive or negative (string)
  #' @param plot: Parameter indicating whether plots should be made (boolean)
  #' @param thresh: Value for noise level threshold (integer)
  #' @param width: Value for width of plot (integer)
  #' @param height: Value for height of plot (integer)
  #'
  #' @return save output to file

  sample_name <- colnames(sample_avgtechrepl)[1]

  # turn dataframe with intensities into a named list
  ints_fullrange <- as.vector(sample_avgtechrepl)
  names(ints_fullrange) <- rownames(sample_avgtechrepl)

  # initialise list to store results for all peaks
  allpeaks_values <- list("mean" = NULL, "area" = NULL, "nr" = NULL,
                          "min" = NULL, "max" = NULL, "qual" = NULL, "spikes" = 0)

  # look for m/z range for all peaks
  allpeaks_values <- search_mzrange(ints_fullrange, allpeaks_values, int_factor, scale, resol,
                                    outdir, sample_name, scanmode,
                                    plot, width, height, thresh)

  # turn the list into a dataframe
  outlist_persample <- NULL
  outlist_persample <- cbind("samplenr" = allpeaks_values$nr,
                             "mzmed.pkt" = allpeaks_values$mean,
                             "fq" = allpeaks_values$qual,
                             "mzmin.pkt" = allpeaks_values$min,
                             "mzmax.pkt" = allpeaks_values$max,
                             "height.pkt" = allpeaks_values$area)

  # remove peaks with height = 0
  outlist_persample <- outlist_persample[outlist_persample[, "height.pkt"] != 0, ]

  # save output to file
  save(outlist_persample, file = paste0(sample_name, "_", scanmode, ".RData"))

  # generate text output to log file on number of spikes for this sample
  # spikes are peaks that are too narrow, e.g. 1 data point
  cat(paste("There were", allpeaks_values$spikes, "spikes"))
}

