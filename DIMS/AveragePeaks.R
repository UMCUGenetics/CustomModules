# load required packages
library("dplyr")
library("argparse")

parser <- ArgumentParser(description = "AveragePeaks")

parser$add_argument("--sample_id", dest = "sample_name",
                    help = "Name of a biological sample", required = TRUE)
parser$add_argument("--tech_reps", dest = "tech_reps",
                    help = "Names of the technical replicates belonging to the biological sample", required = TRUE)
parser$add_argument("--scanmode", dest = "scanmode",
                    help = "Scan mode (either posiive or negative)", required = TRUE)
parser$add_argument("--preprocessing_scripts_dir", dest = "preprocessing_scripts_dir",
                    help = "File path to the directory containing functions used", required = TRUE)

args <- parser$parse_args()

# define parameters
sample_name <- args$sample_id
tech_reps <- strsplit(args$tech_reps, ";")[[1]]

# load in function scripts
source(paste0(args$preprocessing_scripts_dir, "average_peaks_functions.R"))

# Initialize per sample
peaklist_allrepl <- NULL
nr_repl_persample <- 0
averaged_peaks <- matrix(0, nrow = 0, ncol = 6) 
colnames(averaged_peaks) <- c("samplenr", "mzmed.pkt", "fq", "mzmin.pkt", "mzmax.pkt", "height.pkt")

# load RData files of technical replicates belonging to biological sample
for (file_nr in 1:length(tech_reps)) {
  tech_repl_file <- paste0(tech_reps[file_nr], "_", args$scanmode, ".RData")
  tech_repl <- get(load(tech_repl_file))
  # combine data for all technical replicates
  peaklist_allrepl <- rbind(peaklist_allrepl, tech_repl)
}
# sort on mass
peaklist_allrepl_df <- as.data.frame(peaklist_allrepl)
peaklist_allrepl_df$mzmed.pkt <- as.numeric(peaklist_allrepl_df$mzmed.pkt) 
peaklist_allrepl_df$height.pkt <- as.numeric(peaklist_allrepl_df$height.pkt) 
peaklist_allrepl_sorted <- peaklist_allrepl_df %>% arrange(mzmed.pkt)

# average over technical replicates
averaged_peaks <- average_peaks_per_sample(peaklist_allrepl_sorted, sample_name)
save(averaged_peaks, file = paste0("AvgPeaks_", sample_name, "_", args$scanmode, ".RData"))
