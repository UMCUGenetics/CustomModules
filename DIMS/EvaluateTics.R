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
trim_params_filepath <- cmd_args[6]
thresh2remove <- 1000000000

# load init_file: contains repl_pattern
load(init_file)

# load trim_params_file: contains trim_left_neg, trim_left_pos, trim_right_neg & trim_right_pos
load(trim_params_filepath)

# lower the threshold for non Plasma matrices
if (dims_matrix != "Plasma") {
  thresh2remove <- 1000000000
}
# lower the threshold for DBS or high m/z specific
if (dims_matrix == "DBS") {
  thresh2remove <- 500000000
}
if (highest_mz > 700) {
  thresh2remove <- 1000000
}

# find out which technical replicates are below the threshold
remove_tech_reps <- find_bad_replicates(repl_pattern, thresh2remove)
print(remove_tech_reps)

# negative scan mode
remove_neg <- remove_tech_reps$neg
repl_pattern_filtered <- remove_from_repl_pattern(remove_neg, repl_pattern, nr_replicates)
save(repl_pattern_filtered, file = "negative_repl_pattern.RData")

# positive scan mode
remove_pos <- remove_tech_reps$pos
repl_pattern_filtered <- remove_from_repl_pattern(remove_pos, repl_pattern, nr_replicates)
save(repl_pattern_filtered, file = "positive_repl_pattern.RData")

# get an overview of suitable technical replicates for both scan modes
allsamples_techreps_neg <- get_overview_tech_reps(repl_pattern_filtered, "negative")
allsamples_techreps_pos <- get_overview_tech_reps(repl_pattern_filtered, "positive")
allsamples_techreps_both_scanmodes <- rbind(allsamples_techreps_pos, allsamples_techreps_neg)
write.table(allsamples_techreps_both_scanmodes,
            file = "replicates_per_sample.txt",
            col.names = FALSE,
            row.names = FALSE,
            sep = ","
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
    # read file with retention time, intensity and dims_threshold values
    repl1_nr <- read.table(paste0(tech_reps[file_nr], "_TIC.txt"))
    # get threshold values per technical replicate
    dims_thresh_pos <- repl1_nr[1, "threshold"]
    dims_thresh_neg <- repl1_nr[nrow(repl1_nr), "threshold"]
    # for replicates with bad TIC, determine what color the border of the plot should be
    bad_color_pos <- tech_reps[file_nr] %in% remove_pos
    bad_color_neg <- tech_reps[file_nr] %in% remove_neg
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
      geom_vline(xintercept = trim_left_pos, col = "red", linetype = 2, linewidth = 0.3) +
      geom_vline(xintercept = trim_right_pos, col = "red", linetype = 2, linewidth = 0.3) +
      geom_vline(xintercept = trim_left_neg, col = "red", linetype = 2, linewidth = 0.3) +
      geom_vline(xintercept = trim_right_neg, col = "red", linetype = 2, linewidth = 0.3) +
      geom_segment(x = trim_left_pos, y = dims_thresh_pos, xend = trim_right_pos, yend = dims_thresh_pos, colour = "green", lty = 2) + 
      geom_segment(x = trim_left_neg, y = dims_thresh_neg, xend = trim_right_neg, yend = dims_thresh_neg, colour = "blue", lty = 2) + 
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


