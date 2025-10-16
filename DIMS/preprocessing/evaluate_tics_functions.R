# EvaluateTics functions
find_bad_replicates <- function(repl_pattern, thresh2remove) {
  #' Find technical replicates with a total intensity below a threshold
  #'
  #' @param repl_pattern: List of samples with corresponding technical replicates (strings)
  #' @param thresh2remove: Threshold value for acceptance or rejection of total intensity (integer)
  #'
  #' @return remove_tech_reps: Array of rejected technical replicates (strings)
  
  remove_pos <- NULL
  remove_neg <- NULL
  cat("Pklist sum threshold to remove technical replicate:", thresh2remove, "\n")
  for (sample_nr in 1:length(repl_pattern)) {
    tech_reps <- as.vector(unlist(repl_pattern[sample_nr]))
    for (file_nr in 1:length(tech_reps)) {
      load(paste0(tech_reps[file_nr], ".RData"))
      cat("\n\nParsing", tech_reps[file_nr])
      # positive scan mode: determine whether sum of intensities is above threshold
      cat("\n\tPositive peak_list sum", sum(peak_list$pos[, 1]))
      if (sum(peak_list$pos[, 1]) < thresh2remove) {
        cat(" ... Removed")
        remove_pos <- c(remove_pos, tech_reps[file_nr])
      }
      # negative scan mode: determine whether sum of intensities is above threshold
      cat("\n\tNegative peak_list sum", sum(peak_list$neg[, 1]))
      if (sum(peak_list$neg[, 1]) < thresh2remove) {
        cat(" ... Removed")
        remove_neg <- c(remove_neg, tech_reps[file_nr])
      }
    }
  }
  cat("\n")
  # write information on miss_infusions for both scan modes
  write.table(remove_pos,
              file = paste0("miss_infusions_positive.txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t"
  )
  write.table(remove_neg,
              file = paste0("miss_infusions_negative.txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t"
  )
  
  # combine removed technical replicates from pos and neg
  remove_tech_reps <- list(pos = remove_pos, neg = remove_neg)
  return(remove_tech_reps)
}

remove_from_repl_pattern <- function(bad_samples, repl_pattern, nr_replicates) {
  #' Remove technical replicates with insufficient quality from a biological sample
  #'
  #' @param bad_samples: Array of technical replicates of insufficient quality (strings)
  #' @param repl_pattern: List of samples with corresponding technical replicates (strings)
  #' @param nr_replicates: Number of technical replicates per biological sample (integer)
  #'
  #' @return repl_pattern_filtered: list of technical replicates of sufficient quality (strings)
  
  # collect list of samples to remove from replication pattern
  remove_from_group <- NULL
  for (sample_nr in 1:length(repl_pattern)){
    repl_pattern_1sample <- repl_pattern[[sample_nr]]
    remove <- NULL
    for (file_nr in 1:length(repl_pattern_1sample)) {
      if (repl_pattern_1sample[file_nr] %in% bad_samples) {
        remove <- c(remove, file_nr)
      }
    }
    if (length(remove) == nr_replicates) {
      remove_from_group <- c(remove_from_group, sample_nr)
    }
    if (!is.null(remove)) {
      repl_pattern[[sample_nr]] <- repl_pattern[[sample_nr]][-remove]
    }
  }
  if (length(remove_from_group) != 0) {
    repl_pattern_filtered <- repl_pattern[-remove_from_group]
  } else {
    repl_pattern_filtered <- repl_pattern
  }
  return(repl_pattern_filtered)
}

get_overview_tech_reps <- function(repl_pattern_filtered, scanmode) {
  #' Create an overview of technical replicates with sufficient quality from a biological sample
  #'
  #' @param repl_pattern_filtered: List of samples with corresponding technical replicates (strings)
  #' @param scanmode: Scan mode "positive" or "negative" (string)
  #'
  #' @return allsamples_techreps_scanmode: Matrix of technical replicates of sufficient quality (strings)
  
  allsamples_techreps_scanmode <- matrix("", ncol = 3, nrow = length(repl_pattern_filtered))
  for (sample_nr in 1:length(repl_pattern_filtered)) {
    allsamples_techreps_scanmode[sample_nr, 1] <- names(repl_pattern_filtered)[sample_nr]
    allsamples_techreps_scanmode[sample_nr, 2] <- paste0(repl_pattern_filtered[[sample_nr]], collapse = ";")
  }
  allsamples_techreps_scanmode[, 3] <- scanmode
  return(allsamples_techreps_scanmode)
}

