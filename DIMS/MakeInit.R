#!/usr/bin/env Rscript
## adapted from makeInit in old pipeline

# used for when init.RData has to be created manually
# arg1 : path to sampleNames.txt or whatever the name of the samplesheet txt file may be
# arg2 : amount of technical replicates (usually 3)

args <- commandArgs(trailingOnly=TRUE)
sample_sheet <- read.csv(args[1], sep="\t") 
nr_replicates <- as.numeric(args[2])

sampleNames <- trimws(as.vector(unlist(sample_sheet[1])))
nr_sampgrps <- length(sampleNames)/nr_replicates
groupNames <- trimws(as.vector(unlist(sample_sheet[2])))
groupNames <- gsub('[^-.[:alnum:]]', '_', groupNames)
groupNamesUnique <- unique(groupNames)
#groupNamesNotUnique <- groupNames[duplicated(groupNames)]

repl.pattern <- c()
for (sampgrp in 1:nr_sampgrps) {
  tmp <- c()
  for (repl in nr_replicates:1) {
    index <- ((sampgrp*nr_replicates) - repl) + 1
    tmp <- c(tmp, sampleNames[index])
  }
  repl.pattern <- c(repl.pattern, list(tmp))
}

names(repl.pattern) <- groupNamesUnique

# just to preview
head(repl.pattern)

save(repl.pattern, file="init.RData")
