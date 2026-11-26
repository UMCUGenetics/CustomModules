# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)

db_file <- cmd_args[1]
breaks_file <- cmd_args[2]

load(db_file)
load(breaks_file)

# get minimum and maximum m/z in dataset
min_mz <- round(breaks_fwhm[1])
max_mz <- round(breaks_fwhm[length(breaks_fwhm)])

# Select HMDB plus adducts and isotopes for each scan mode 
scanmodes <- c("positive", "negative")
for (scanmode in scanmodes) {
  if (scanmode == "negative") {
    column_label <- "MNeg"
    HMDB_add_iso <- HMDB_add_iso.Neg
  } else if (scanmode == "positive") {
    column_label <- "Mpos"
    HMDB_add_iso <- HMDB_add_iso.Pos
  }

  # filter on mass range in dataset
  HMDB_mzrange <- HMDB_add_iso[(HMDB_add_iso[, column_label] >= min_mz & HMDB_add_iso[, column_label] <= max_mz), ]
  # remove adducts and isotopes, put internal standard at the beginning
  HMDB_add <- HMDB_mzrange[grep("HMDB", rownames(HMDB_mzrange), fixed = TRUE), ]
  HMDB_main <- HMDB_add[-grep("_", rownames(HMDB_add), fixed = TRUE), ]
  # sort on m/z value
  HMDB_main <- HMDB_main[order(HMDB_main[, column_label]), ]
  
  # generate hmdb parts of 1000 lines each
  nr_parts <- ceiling(nrow(HMDB_main) / 1000)
  start_index <- 1
  for (part_index in 1:nr_parts) {
    end_index <- min((start_index + 999), nrow(HMDB_main))
    outlist_part <- HMDB_main[start_index:end_index, ]
    save(outlist_part, file = paste0(scanmode, "_hmdb_main.", part_index, ".RData"))
    start_index = start_index + 1000
  }
}
