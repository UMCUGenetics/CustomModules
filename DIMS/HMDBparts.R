# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)

db_file <- cmd_args[1] 
breaks_file <- cmd_args[2] 
standard_run <- cmd_args[3] 

# load file with binning breaks 
load(breaks_file)
min_mz <- round(breaks_fwhm[1])
max_mz <- round(breaks_fwhm[length(breaks_fwhm)])

# In case of a standard run use external HMDB parts
# m/z is approximately 70 to 600: set limits between 68-71 for min and 599-610 for max
if (standard_run == "yes" & min_mz > 68 & min_mz < 71 & max_mz > 599 & max_mz < 610) {
  # skip generating HMDB parts; copy them from hmdb_parts_path
  hmdb_parts_path <- cmd_args[4] 
  # find all files containing hmdb in file name
  hmdb_parts <- list.files(hmdb_parts_path, pattern = "hmdb") 
  for (hmdb_file in hmdb_parts) {
    file.copy(paste(hmdb_parts_path, hmdb_file, sep = "/"), "./", recursive = TRUE)
  }
} else { 
  # generate HMDB parts in case of non-standard mz range
  load(db_file)

  # determine segments of m/z for HMDB parts; smaller parts for m/z < 100
  mz_segments <- c()
  segment_start <- min_mz
  segment_end <- min_mz + 5
  while (segment_end < max_mz) {
    if (segment_start < 100) {
      mz_segments <- c(mz_segments, segment_start)
      segment_start <- segment_start + 5
      segment_end <- segment_end + 5
    } else {
      mz_segments <- c(mz_segments, segment_start)
      segment_start <- segment_start + 10
      segment_end <- segment_end + 10
    }
  }
  #last segment
  mz_segments <- c(mz_segments, max_mz)
  
  scanmodes <- c("positive", "negative")
  for (scanmode in scanmodes) {
    if (scanmode == "negative") {
      column_label <- "MNeg"
      hmdb_add_iso <- HMDB_add_iso.Neg
    } else if (scanmode == "positive") {
      column_label <- "Mpos"
      hmdb_add_iso <- HMDB_add_iso.Pos
    }

    # filter mass range meassured
    hmdb_add_iso = hmdb_add_iso[which(hmdb_add_iso[ , column_label] >= min_mz &
                                      hmdb_add_iso[ , column_label] <= max_mz), ]

    # sort on mass
    sorted_hmdb_add_iso <- hmdb_add_iso[order(as.numeric(hmdb_add_iso[ , column_label])),]
    nr_rows <- dim(sorted_hmdb_add_iso)[1]

    # create parts and save to file
    for (mz_part_index in 1:(length(mz_segments) - 1)) {
      mz_start <- mz_segments[mz_part_index]
      mz_end <- mz_segments[mz_part_index + 1]
      outlist_part <- sorted_hmdb_add_iso[sorted_hmdb_add_iso[ , column_label] > mz_start &
                                          sorted_hmdb_add_iso[ , column_label] <= mz_end, ]
      save(outlist_part, file = paste0(scanmode, "_hmdb.", mz_part_index, ".RData"))
    }
  }
}
