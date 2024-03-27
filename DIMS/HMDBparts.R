#!/usr/bin/Rscript
# adapted from hmdb_parts.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)

db_file <- cmd_args[1] 
breaks_file <- cmd_args[2] 
standard_run <- cmd_args[4] 

# load file with binning breaks 
load(breaks_file)
min_mz <- round(breaks_fwhm[1])
max_mz <- round(breaks_fwhm[length(breaks_fwhm)])

# In case of a standard run (m/z 69-606) use external HMDB parts
if (standard_run == "yes" & min_mz > 68 & min_mz < 71 & max_mz > 599 & max_mz < 610) {
  # skip generating HMDB parts
  hmdb_parts_path <- cmd_args[3] 
  # find all files containing hmdb in file name
  hmdb_parts <- list.files(hmdb_parts_path, pattern = "hmdb") 
  for (hmdb_file in hmdb_parts) {
    file.copy(paste(hmdb_parts_path, hmdb_file, sep = "/"), "./", recursive = TRUE)
  }
} else { 
  # generate HMDB parts in case of non-standard mz range
  load(db_file)
  ppm <- as.numeric(cmd_args[5])

  scanmodes <- c("positive", "negative")
  for (scanmode in scanmodes) {
    if (scanmode == "negative") {
      column_label <- "MNeg"
      HMDB_add_iso <- HMDB_add_iso.Neg
    } else if (scanmode == "positive") {
      column_label <- "Mpos"
      HMDB_add_iso <- HMDB_add_iso.Pos
    }

    # filter mass range meassured
    HMDB_add_iso = HMDB_add_iso[which(HMDB_add_iso[ , column_label] >= breaks_fwhm[1] &
                                      HMDB_add_iso[ , column_label] <= breaks_fwhm[length(breaks_fwhm)]), ]

    # sort on mass
    outlist <- HMDB_add_iso[order(as.numeric(HMDB_add_iso[ , column_label])),]
    nr_rows <- dim(outlist)[1]

    # maximum number of rows per file
    sub <- 20000 
    end <- 0
    last_line <- sub
    check <- 0
    outlist_part <- NULL

    # create parts and save to file
    if (nr_rows < sub) {
      outlist_part <- outlist
      save(outlist_part, file = paste0(scanmode, "_hmdb.1.RData"))
    } else if (nr_rows >= sub & (floor(nr_rows / sub) - 1) >= 2) {
      for (i in 2:floor(nr_rows / sub) - 1) {
        start <- -(sub - 1) + i * sub
        end <- i * sub

        if (i > 1){
          outlist_i = outlist[c(start:end),]
          nr_moved = 0
          # Use ppm to replace border to avoid cut within peak group
          while ((as.numeric(outlist_i[1, column_label]) - as.numeric(outlist_part[last_line, column_label])) * 1e+06 / 
                  as.numeric(outlist_i[1, column_label]) < ppm) {
            outlist_part <- rbind(outlist_part, outlist_i[1, ])
            outlist_i <- outlist_i[-1, ]
            nr_moved <- nr_moved + 1
          }

          save(outlist_part, file = paste(scanmode, "_", paste("hmdb", i-1, "RData", sep = "."), sep = ""))
          check <- check + dim(outlist_part)[1]

          outlist_part <- outlist_i
          last_line <- dim(outlist_part)[1]

        } else {
          outlist_part <- outlist[c(start:end),]
        }
      }

      start <- end + 1
      end <- nr_rows
      outlist_i <- outlist[c(start:end), ]
      nr_moved <- 0

      if (!is.null(outlist_part)) {
        # Calculate ppm and replace border, avoid cut within peak group
        while ((as.numeric(outlist_i[1, column_label]) - as.numeric(outlist_part[last_line, column_label])) * 1e+06 / 
                as.numeric(outlist_i[1, column_label]) < ppm) {
          outlist_part <- rbind(outlist_part, outlist_i[1, ])
          outlist_i <- outlist_i[-1, ]
          nr_moved <- nr_moved + 1
        }

        save(outlist_part, file = paste0(scanmode, "_hmdb_", i, ".RData"))
        check <- check + dim(outlist_part)[1]
      }

      outlist_part <- outlist_i
      save(outlist_part, file = paste0(scanmode, "_hmdb_", i + 1, ".RData"))
      check <- check + dim(outlist_part)[1]
    }
  }
}

