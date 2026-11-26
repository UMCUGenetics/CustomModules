## adapted from checkOverlap.R
check_overlap <- function(range1, range2) {
  #' Modify range1 and range2 in case of overlap
  #'
  #' @param range1: Vector of m/z values for first peak (float)
  #' @param range2: Vector of m/z values for second peak (float)
  #'
  #' @return new_ranges: list of two ranges (list)
  
  # Check for overlap
  if (length(intersect(range1, range2)) == 2) {
    if (length(range1) >= length(range2)) {
      range1 <- range1[-length(range1)]  
    } else {
      range2 <- range2[-1]
    }
  } else if (length(intersect(range1, range2)) == 3) {
    range1 <- range1[-length(range1)]  
    range2 <- range2[-1]
  }
  new_ranges <- list("range1" = range1, "range2" = range2)
  return(new_ranges)
}

