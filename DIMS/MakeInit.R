## adapted from makeInit in old pipeline

# create object for run with constructor function
Run <- function(run_name, date, samples, nr_replicates) {
  run_obj <- list(
    run_name = run_name,
    date = date,
    samples = samples,
    nr_replicates = nr_replicates
  )
  class(run_obj) <- "Run"
  return(run_obj)
}

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

TechnicalReplicate <- function(sample_name, peaklist_pos, peaklist_neg, raw_file) {
  replicate_obj <- list(
    sample_name = sample_name,
    peaklist_pos = NULL,
    peaklist_neg = NULL,
    raw_file = raw_file
  )
  class(replicate_obj) <- "TechnicalReplicate"
  return(replicate_obj)
}

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

# make a list of samples
list_samples <- lapply(1:nr_samples, function(sample_nr) {
  indices <- ((sample_nr - 1) * nr_replicates + 1):(sample_nr * nr_replicates)
  sample_name <- sample_names_unique[sample_nr]
  
  # make a list of technical replicates
  list_tech_reps <- lapply(indices, function(index) {
    TechnicalReplicate(sample_name = sample_name, 
                       peaklist_pos = NULL, 
                       peaklist_neg = NULL, 
                       raw_file = file_names[index])
  })
  names(list_tech_reps) <- file_names[indices]
  
  # make a Sample object
  Sample(sample_name = sample_name, 
         technical_replicates = list_tech_reps, 
         averaged_peaklist_pos = NULL, 
         averaged_peaklist_neg = NULL)
})

names(list_samples) <- sample_names_unique

current_run <- Run(run_name = run_name, samples = list_samples, nr_replicates = nr_replicates, date = Sys.Date())
save(current_run, file = "run_info.RData")


## Old code, unnecessary with OOP ##

# generate the replication pattern
repl_pattern <- lapply(1:nr_samples, function(sample) {
  indices <- ((sample - 1) * nr_replicates + 1):(sample * nr_replicates)
  file_names[indices]
})
names(repl_pattern) <- sample_names_unique

# preview the replication pattern
print(tail(repl_pattern))

save(repl_pattern, file = "init.RData")
