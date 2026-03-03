# load required packages
library("argparse")

parser <- ArgumentParser(description = "MakeInit")

parser$add_argument("--samplesheet", dest = "sample_sheet",
                    help = "Samplesheet txt file", required = TRUE)
parser$add_argument("--nr_replicates", dest = "nr_replicates", type = "integer",
                    help = "Number of replicates, numeric value", required = TRUE)

args <- parser$parse_args()

sample_sheet <- read.csv(args$sample_sheet, sep = "\t")
sample_names <- trimws(as.vector(unlist(sample_sheet[1])))
nr_replicates <- args$nr_replicates
nr_sample_groups <- length(sample_names) / nr_replicates
group_names <- trimws(as.vector(unlist(sample_sheet[2])))
group_names <- gsub("[^-.[:alnum:]]", "_", group_names)
group_names_unique <- unique(group_names)

# generate the replication pattern
repl_pattern <- c()
for (sample_group in 1:nr_sample_groups) {
  replicates_persample <- c()
  for (repl in nr_replicates:1) {
    index <- ((sample_group * nr_replicates) - repl) + 1
    replicates_persample <- c(replicates_persample, sample_names[index])
  }
  repl_pattern <- c(repl_pattern, list(replicates_persample))
}

names(repl_pattern) <- group_names_unique

# preview the replication pattern
print(tail(repl_pattern))

save(repl_pattern, file = "init.RData")
