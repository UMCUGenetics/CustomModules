## adapted from makeInit in old pipeline

# define parameters
args <- commandArgs(trailingOnly = TRUE)

sample_sheet <- read.csv(args[1], sep = "\t")
nr_replicates <- as.numeric(args[2])
run_name <- args[3]

file_names <- base::trimws(as.vector(unlist(sample_sheet[1])))
nr_samples <- length(file_names) / nr_replicates
sample_names <- base::trimws(as.vector(unlist(sample_sheet[2])))
sample_names <- gsub("[^-.[:alnum:]]", "_", sample_names)
sample_names_unique <- unique(sample_names)

# generate the replication pattern
repl_pattern <- lapply(1:nr_samples, function(sample) {
  indices <- ((sample -1) * nr_replicates+ 1):(sample * nr_replicates)
  file_names[indices]
})
names(repl_pattern) <- sample_names_unique

# preview the replication pattern
print(tail(repl_pattern))

# create object for run with constructor function
Run <- function(run_name, date, sample_names, nr_replicates) {
  run_obj <- list(
    run_name = run_name,
    date = date,
    sample_names = sample_names,
    nr_replicates = nr_replicates
  )
  class(run_obj) <- "Run"
  return(run_obj)
}

current_run <- Run(run_name = run_name, samples = group_names, nr_replicates = nr_replicates)
save(current_run, file = "run_info.RData")

Sample <- function(sample_name, technical_replicates, averaged_peaklist_pos, averaged_peaklist_neg) {
  sample_obj <- list(
    sample_name = sample_name,
    technical_replicates = technical_replicates,
    averaged_peaklist_pos = NULL,
    averaged_peaklist_neg = NULL
  )
  class(sample_obj) <- "Sample"
  return(sample_obj)
}

TechnicalReplicate <- function(sammple_name, peaklist_pos, peaklist_neg) {
  replicate_obj <- list(
    sample_name = sample_name,
    peaklist_pos = NULL,
    peaklist_neg = NULL
  )
  class(replicate_obj) <- "TechnicalReplicate"
  return(replicate_obj)
}


save(repl_pattern, file = "init.RData")
