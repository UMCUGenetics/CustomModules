#!/usr/bin/Rscript
## adapted from 11-runSumAdducts.R

# define parameters
cmd_args <- commandArgs(trailingOnly = TRUE)

hmdbpart_main_file <- cmd_args[1]
scripts_dir <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])

sum_adducts <- function(peaklist, theor_mz, grpnames_long, adducts, batch_number, scanmode, outdir, z_score) {
  hmdb_codes <- rownames(theor_mz)
  hmdb_names <- theor_mz[, 1, drop = FALSE]
  hmdb_names[] <- lapply(hmdb_names, as.character)

  # remove isotopes
  index <- grep("HMDB", hmdb_codes, fixed = TRUE)
  hmdb_codes <- hmdb_codes[index]
  hmdb_names <- hmdb_names[index, ]
  index <- grep("_", rownames(hmdb_codes), fixed = TRUE)
  if (length(index) > 0) hmdb_codes <- hmdb_codes[-index]
  if (length(index) > 0) hmdb_names <- hmdb_names[-index]

  # negative
  names <- NULL
  adductsum <- NULL
  names_long <- NULL

  if (length(hmdb_codes) != 0) {
    for (i in 1:length(hmdb_codes)) {
      compound <- hmdb_codes[i]
      compound_plus <- c(compound, paste(compound, adducts, sep = "_"))

      metab <- unlist(lapply(peaklist$HMDB_code, function(x) {
        (length(intersect(unlist(strsplit(as.vector(x), ";")), compound_plus)) > 0)
      }))

      total <- c()

      # get the intensities for selected metabolite.
      if (z_score == 1) {
        int_cols_C <- grep("C", colnames(peaklist)[1:which(colnames(peaklist) == "avg.ctrls")])
        int_cols_P <- grep("P", colnames(peaklist)[1:which(colnames(peaklist) == "avg.ctrls")])
        int_cols <- c(int_cols_C, int_cols_P)
        ints <- peaklist[metab, int_cols]
      } else {
        ints <- peaklist[metab, c(7:(length(grpnames_long) + 6))]
      }
      total <- apply(ints, 2, sum)

      if (sum(total) != 0) {
        names <- c(names, compound)
        adductsum <- rbind(adductsum, total)
        names_long <- c(names_long, hmdb_names[i])
      }
    }

    if (!is.null(adductsum)) {
      rownames(adductsum) <- names
      adductsum <- cbind(adductsum, "HMDB_name" = names_long)
      save(adductsum, file = paste(scanmode, "_", batch_number, "_SummedAdducts.RData", sep = ""))
    }
  }
}

if (grepl("positive_hmdb", hmdbpart_main_file)) {
  scanmode <- "positive"
  # for the adduct sum: include adducts M+Na (1) and M+K (2)
  adducts <- c(1, 2)
} else if (grepl("negative_hmdb", hmdbpart_main_file)) {
  scanmode <- "negative"
  # for the adduct sum: include adduct M+Cl (1)
  adducts <- c(1)
}

# load input files
collect_file <- paste0("outlist_identified_", scanmode, ".RData")
load(collect_file)
repl_file <- paste0(scanmode, "_repl_pattern.RData")
load(repl_file)
outlist_part <- get(load(hmdbpart_main_file))

# get the number from the file name
batch_number <- strsplit(basename(hmdbpart_main_file), ".", fixed = TRUE)[[1]][2]

outlist_total <- unique(outlist_ident)

sum_adducts(outlist_total, outlist_part, names(repl_pattern_filtered), adducts, batch_number, scanmode, outdir, z_score)
