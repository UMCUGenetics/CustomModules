#!/usr/bin/Rscript
# adapted from 3-AverageTechReplicates.R

# load packages
library("ggplot2")
library("gridExtra")

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

init_file <- cmd_args[1]
nr_replicates <- as.numeric(cmd_args[2])
run_name <- cmd_args[3]
dims_matrix <- cmd_args[4]
highest_mz_file <- cmd_args[5]
highest_mz <- get(load(highest_mz_file))
thresh2remove <- 1000000000
dims_thresh <- 100

remove_from_repl_pattern <- function(bad_samples, repl_pattern, nr_replicates) {
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
    repl_pattern <- repl_pattern[-remove_from_group]
  }
  return(list("pattern" = repl_pattern))
}

# get repl_pattern
load(init_file)

# lower the threshold below which a sample will be removed for DBS and for high m/z
if (dims_matrix == "DBS") {
  thresh2remove <- 500000000
}
if (highest_mz > 700) {
  thresh2remove <- 1000000
}

# remove technical replicates which are below the threshold
remove_neg <- NULL
remove_pos <- NULL
cat("Pklist sum threshold to remove technical replicate:", thresh2remove, "\n")
for (sample_nr in 1:length(repl_pattern)) {
  tech_reps <- as.vector(unlist(repl_pattern[sample_nr]))
  tech_reps_array_pos <- NULL
  tech_reps_array_neg <- NULL
  sum_neg <- 0
  sum_pos <- 0
  nr_pos <- 0
  nr_neg <- 0
  for (file_nr in 1:length(tech_reps)) {
    load(paste(tech_reps[file_nr], ".RData", sep = ""))
    cat("\n\nParsing", tech_reps[file_nr])
    # negative scanmode
    cat("\n\tNegative peak_list sum", sum(peak_list$neg[, 1]))
    if (sum(peak_list$neg[, 1]) < thresh2remove) {
      cat(" ... Removed")
      remove_neg <- c(remove_neg, tech_reps[file_nr])
    } else {
      nr_neg <- nr_neg + 1
      sum_neg <- sum_neg + peak_list$neg
    }
    tech_reps_array_neg <- cbind(tech_reps_array_neg, peak_list$neg)
    # positive scanmode
    cat("\n\tPositive peak_list sum", sum(peak_list$pos[, 1]))
    if (sum(peak_list$pos[, 1]) < thresh2remove) {
      cat(" ... Removed")
      remove_pos <- c(remove_pos, tech_reps[file_nr])
    } else {
      nr_pos <- nr_pos + 1
      sum_pos <- sum_pos + peak_list$pos
    }
    tech_reps_array_pos <- cbind(tech_reps_array_pos, peak_list$pos)
  }
  # save to file  
  if (nr_neg != 0) {
    sum_neg[, 1] <- sum_neg[, 1] / nr_neg
    colnames(sum_neg) <- names(repl_pattern)[sample_nr]
    save(sum_neg, file = paste0(names(repl_pattern)[sample_nr], "_neg_avg.RData"))
  }
  if (nr_pos != 0) {
    sum_pos[, 1] <- sum_pos[, 1] / nr_pos
    colnames(sum_pos) <- names(repl_pattern)[sample_nr]
    save(sum_pos, file = paste0(names(repl_pattern)[sample_nr], "_pos_avg.RData"))
  }
}

pattern_list <- remove_from_repl_pattern(remove_neg, repl_pattern, nr_replicates)
repl_pattern_filtered <- pattern_list$pattern
save(repl_pattern_filtered, file = "negative_repl_pattern.RData")
write.table(
  remove_neg, 
  file = "miss_infusions_negative.txt", 
  row.names = FALSE, 
  col.names = FALSE, 
  sep = "\t"
)

pattern_list <- remove_from_repl_pattern(remove_pos, repl_pattern, nr_replicates)
repl_pattern_filtered <- pattern_list$pattern
save(repl_pattern_filtered, file = "positive_repl_pattern.RData")
write.table(
  remove_pos, 
  file = "miss_infusions_positive.txt", 
  row.names = FALSE, 
  col.names = FALSE, 
  sep = "\t"
)

## generate TIC plots
# get all txt files
tic_files <- list.files("./", full.names = TRUE, pattern = "*TIC.txt")
all_samps <- sub("_TIC\\..*$", "", basename(tic_files))

# determine maximum intensity
highest_tic_max <- 0
for (file in tic_files) {
  tic <- read.table(file)
  this_tic_max <- max(tic$tic_intensity)
  if (this_tic_max > highest_tic_max) {
    highest_tic_max <- this_tic_max
    max_sample <- sub("_TIC\\..*$", "", basename(file))
  }
}

# create a list with information for all TIC plots
tic_plot_list <- list()
plot_nr <-  0
for (sample_nr in c(1:length(repl_pattern))) {
  tech_reps <- as.vector(unlist(repl_pattern[sample_nr]))
  sample_name <- names(repl_pattern)[sample_nr]
  for (file_nr in 1:length(tech_reps)) {
    plot_nr <- plot_nr + 1
    # repl1_nr <- read.table(paste(paste(outdir, "2-pklist/", sep = "/"), tech_reps[file_nr], "_TIC.txt", sep = ""))
    repl1_nr <- read.table(paste0(tech_reps[file_nr], "_TIC.txt"))
    bad_color_pos <- tech_reps[file_nr] %in% remove_pos[[1]]
    bad_color_neg <- tech_reps[file_nr] %in% remove_neg[[1]]
    if (bad_color_neg & bad_color_pos) {
      plot_color <- "#F8766D"
    } else if (bad_color_pos) {
      plot_color <- "#ED8141"
    } else if (bad_color_neg) {
      plot_color <- "#BF80FF"
    } else {
      plot_color <- "white"
    }
    tic_plot <- ggplot(repl1_nr, aes(retention_time, tic_intensity)) +
      geom_line(linewidth = 0.3) +
      geom_hline(yintercept = highest_tic_max, col = "grey", linetype = 2, linewidth = 0.3) +
      labs(x = "t (s)", y = "tic_intensity", title = paste0(tech_reps[file_nr], "  ||  ", sample_name)) +
      theme(plot.background = element_rect(fill = plot_color),
            axis.text = element_text(size = 4),
            axis.title = element_text(size = 4),
            plot.title = element_text(size = 6))
    tic_plot_list[[plot_nr]] <- tic_plot
  }
}

# create a layout matrix dependent on number of replicates
layout <- matrix(1:(10 * nr_replicates), 10, nr_replicates, TRUE)
# put TIC plots in matrix
tic_plot_pdf <- marrangeGrob(
  grobs = tic_plot_list,
  nrow = 10, ncol = nr_replicates,
  layout_matrix = layout,
  top = quote(paste(
    "TICs of run", run_name,
    " \n colors: red = both modes misinjection, orange = pos mode misinjection, purple = neg mode misinjection \n ",
    g, "/", npages
  ))
)

# save to file
ggsave(filename = paste0(run_name, "_TICplots.pdf"),
       tic_plot_pdf, width = 21, height = 29.7, units = "cm")
