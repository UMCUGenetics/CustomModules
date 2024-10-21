## adapted from makeInit in old pipeline

# define parameters
args <- commandArgs(trailingOnly = TRUE)

#sample_sheet <- read.csv(args[1], sep = "\t")
sample_sheet <- read.csv(args[1], sep = "\t", nrows = 8)  # For debugging purposes
nr_replicates <- as.numeric(args[2])

sample_names <- trimws(as.vector(unlist(sample_sheet[1])))
nr_sample_groups <- length(sample_names) / nr_replicates
group_names <- trimws(as.vector(unlist(sample_sheet[2])))
group_names <- gsub("[^-.[:alnum:]]", "_", group_names)
group_names_unique <- unique(group_names)

# generate the replication pattern
repl_pattern <- c()
for (sample_group in 1:nr_sample_groups) {
  tmp <- c()
  for (repl in nr_replicates:1) {
    index <- ((sample_group * nr_replicates) - repl) + 1
    tmp <- c(tmp, sample_names[index])
  }
  repl_pattern <- c(repl_pattern, list(tmp))
}

names(repl_pattern) <- group_names_unique

# preview the replication pattern
print(tail(repl_pattern))

save(repl_pattern, file = "init.RData")
