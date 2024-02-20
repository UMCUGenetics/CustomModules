#!/usr/bin/env Rscript
## adapted from makeInit in old pipeline

args          <- commandArgs(trailingOnly = TRUE)
sample_sheet  <- read.csv(args[1], sep = "\t")
nr_replicates <- as.numeric(args[2])

sample_names  <- trimws(as.vector(unlist(sample_sheet[1])))
nr_sampgrps   <- length(sampleNames) / nr_replicates
group_names    <- trimws(as.vector(unlist(sample_sheet[2])))
group_names    <- gsub("[^-.[:alnum:]]", "_", group_names)
group_names_unique <- unique(group_names)

repl_pattern <- c()
for (sampgrp in 1:nr_sampgrps) {
  tmp <- c()
  for (repl in nr_replicates:1) {
    index <- ((sampgrp * nr_replicates) - repl) + 1
    tmp <- c(tmp, sampleNames[index])
  }
  repl_pattern <- c(repl_pattern, list(tmp))
}

names(repl_pattern) <- group_names_unique

# preview the replication pattern
print(tail(repl_pattern))

save(repl_pattern, file = "./init.RData")
