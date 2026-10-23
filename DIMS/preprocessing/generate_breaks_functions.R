# GenerateBreaks functions
# functions: get_trim_parameters, get_breaks_for_bins

#' Determine the scans per scanmode which are trimmed off; save trim parameters to file
#'
#' @param scantimes: Scan times in seconds (vector of float)
#' @param polarities: Polarities, positive or negative (vector of string)
#' @param trim: Fraction of scans which are to be discarded (float)
get_trim_parameters <- function(scantimes, polarities, trim = 0.1) {
  # Get time values for positive and negative scans
  pos_times <- scantimes[polarities == "positive"]
  neg_times <- scantimes[polarities == "negative"]

  # Remove scans at the start and end for positive
  trim_left_pos  <- round(pos_times[length(pos_times) * (trim * 1.5)])
  trim_right_pos <- round(pos_times[length(pos_times) * (1 - (trim * 0.5))])

  # Remove scans at the start and end for negative
  trim_left_neg  <- round(neg_times[length(neg_times) * trim])
  trim_right_neg <- round(neg_times[length(neg_times) * (1 - trim)])

  # Save trim parameters to file
  save(trim_left_pos, trim_right_pos, trim_left_neg, trim_right_neg, file = "trim_params.RData")
}

#' create a vector with the breaks in m/z of bins for intensities
#'
#' @param mzrange: Minimum and maximum m/z values (vector of integers)
#' @param resol: Resolution (integer)
get_breaks_for_bins <- function(mzrange, resol = 140000) {
  # Initialize
  breaks_fwhm <- NULL
  breaks_fwhm_avg <- NULL

  # Determine number of segments used to create bins
  nr_segments <- 2 * (mzrange[2] - mzrange[1])
  segments <- seq(from = mzrange[1], to = mzrange[2], length.out = nr_segments + 1)

  # Determine start and end of each bin. fwhm (width of peaks) is assumed to be constant within a segment
  for (segment_index in 1:nr_segments) {
    start_segment <- segments[segment_index]
    end_segment <- segments[segment_index + 1]
    # Determine resolution at given m/z value
    resol_mz <- resol * (1 / sqrt(2) ^ (log2(start_segment / 200)))
    # Determine fwhm (full width at half maximum) of the peaks in this segment
    fwhm_segment <- start_segment / resol_mz
    # Determine the breaks within this segment
    breaks_segment <- seq(from = (start_segment + fwhm_segment), to = end_segment, by = 0.2 * fwhm_segment)
    # Add breaks for this segment to vector with all breaks
    breaks_fwhm <- c(breaks_fwhm, seq(from = (start_segment + fwhm_segment), to = end_segment, by = 0.2 * fwhm_segment))
    # Get a vector of average m/z instead of start value
    delta_mz <- breaks_segment[2] - breaks_segment[1]
    avg_breaks_segment <- breaks_segment + 0.5 * delta_mz
    breaks_fwhm_avg <- c(breaks_fwhm_avg, avg_breaks_segment)
  }

  # Save breaks to file
  save(breaks_fwhm, breaks_fwhm_avg, file = "breaks.fwhm.RData")
}
