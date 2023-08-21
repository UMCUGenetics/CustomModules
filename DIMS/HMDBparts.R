#!/usr/bin/Rscript
# adapted from hmdb_parts.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n")

# outdir <- cmd_args[1]
# scanmode <- cmd_args[2]
db_path <- cmd_args[1] # location of HMDB db file
breaks_filepath <- cmd_args[2] # location of breaks.fwhm.RData
standard_run  <- cmd_args[3] # "yes"

# Cut up entire HMDB into small parts based on the new binning/breaks 

load(breaks_filepath)

# In case of a standard run (m/z 69-606) use external HMDB parts
min_mz <- round(breaks_fwhm[1])
max_mz <- round(breaks_fwhm[length(breaks_fwhm)])
# test if standard mz range is used
if (standard_run == "yes" & min_mz > 68 & min_mz < 71 & max_mz > 600 & max_mz < 610) {
  # skip generating HMDB parts
  use_external_HMDB <- TRUE
  save(use_external_HMDB, file="./using_external_HMDB_parts.RData")
  hmdb_parts_dir <- cmd_args[6] # "/hpc/dbg_mz/production/DIMS/hmdb_preparts/"
  # hmdb_parts <- list.files(hmdb_parts_dir, pattern=scanmode) # only positive or negative files
  hmdb_parts <- list.files(hmdb_parts_dir, pattern=hmdb) # all files containing hmdb in file name
  # save(HMDBstukken, file=paste(outdir, "HMDBstukken.RData", sep="/"))
  for (hmdb_file in hmdb_parts) {
    # file.copy(paste(hmdb_parts_dir, hmdb_file, sep="/"), outdir_hmdb, recursive = TRUE)
    file.copy(paste(hmdb_parts_dir, hmdb_file, sep="/"), "./", recursive = TRUE)
  }
} else { 
  # generate HMDB parts in case of non-standard mz range
  use_external_HMDB <- FALSE
  save(use_external_HMDB, file="./not_using_external_HMDB_parts.RData")
  # db <- cmd_args[3]
  load(db_path)
  ppm <- as.numeric(cmd_args[4])

  scanmodes <- c("positive", "negative")

  for (scanmode in scanmodes) {
    if (scanmode == "negative") {
      column_label <- "MNeg"
      HMDB_add_iso <- HMDB_add_iso.Neg
    } else if (scanmode == "positive") {
      column_label <- "Mpos"
      HMDB_add_iso <- HMDB_add_iso.Pos
    }

    # filter mass range meassured!!!
    HMDB_add_iso = HMDB_add_iso[which(HMDB_add_iso[ ,column_label] >= breaks_fwhm[1] &
                                      HMDB_add_iso[ ,column_label] <= breaks_fwhm[length(breaks_fwhm)]), ]

    # sort on mass
    outlist = HMDB_add_iso[order(as.numeric(HMDB_add_iso[,column_label])),]

    n=dim(outlist)[1]
    sub=20000 # max rows per file
    end=0
    min_1_last=sub
    check=0
    outlist_part=NULL


    if (n < sub) {
      outlist_part <- outlist
      save(outlist_part, file = paste0("./", scanmode, "_hmdb.1.RData"))
    } else {

        if (n >= sub & (floor(n/sub) - 1) >= 2){
          for (i in 2:floor(n/sub) - 1){
            start <- -(sub - 1) + i*sub
            end <- i*sub

            if (i > 1){
              outlist_i = outlist[c(start:end),]

              n_moved = 0

              # Calculate 3ppm and replace border, avoid cut within peakgroup!
              while ((as.numeric(outlist_i[1,column_label]) - as.numeric(outlist_part[min_1_last,column_label]))*1e+06/as.numeric(outlist_i[1,column_label]) < ppm) {
                outlist_part <- rbind(outlist_part, outlist_i[1,])
                outlist_i <- outlist_i[-1,]
                n_moved <- n_moved + 1
              }

              # message(paste("Process", i-1,":", dim(outlist_part)[1]))
              save(outlist_part, file = paste("./", scanmode, "_", paste("hmdb", i-1, "RData", sep="."), sep=""))
              check <- check + dim(outlist_part)[1]

              outlist_part <- outlist_i
              min_1_last <- dim(outlist_part)[1]

            } else {
              outlist_part <- outlist[c(start:end),]
            }
          }
        }

        start <- end + 1
        end <- n
        outlist_i <- outlist[c(start:end),]
        n_moved <- 0

        if (!is.null(outlist_part)) {
          # Calculate 3ppm and replace border, avoid cut within peakgroup!
          while ((as.numeric(outlist_i[1,column_label]) - as.numeric(outlist_part[min_1_last,column_label]))*1e+06/as.numeric(outlist_i[1,column_label]) < ppm) {
            outlist_part = rbind(outlist_part, outlist_i[1,])
            outlist_i = outlist_i[-1,]
            n_moved = n_moved + 1
          }

          # message(paste("Process", i+1-1,":", dim(outlist_part)[1]))
          save(outlist_part, file = paste("./", scanmode, "_", paste("hmdb", i, "RData", sep = "."), sep = ""))
          check <- check + dim(outlist_part)[1]
        }

        outlist_part <- outlist_i
        # message(paste("Process", i+2-1,":", dim(outlist_part)[1]))
        save(outlist_part, file = paste("./", scanmode, "_", paste("hmdb", i + 1, "RData", sep="."), sep=""))
        check <- check + dim(outlist_part)[1]
        cat("\n", "Check", check == dim(outlist)[1])

      }
    }

}

