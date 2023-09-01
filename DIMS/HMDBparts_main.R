#!/usr/bin/Rscript
## adapted from hmdb_part_adductSums.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n")

db_path <- cmd_args[1] # location of HMDB db file
breaks_filepath <- cmd_args[2] # location of breaks.fwhm.RData

load(db_path)
load(breaks_filepath)

# Cut up HMDB minus adducts minus isotopes into small parts 

# load(paste(outdir, "breaks.fwhm.RData", sep = "/"))
# outdir <- paste(outdir, "hmdb_part_adductSums", sep = "/")
# dir.create(outdir, showWarnings = FALSE)

scanmodes <- c("positive", "negative")

for (scanmode in scanmodes) {
  if (scanmode == "negative") {
    column_label <- "MNeg"
    HMDB_add_iso <- HMDB_add_iso.Neg
  } else if (scanmode == "positive") {
    column_label <- "Mpos"
    HMDB_add_iso <- HMDB_add_iso.Pos
  }

  # filter mass range meassured NB: remove the last comma?!
  outlist <- HMDB_add_iso[which(HMDB_add_iso[ ,column_label] >= breaks_fwhm[1] & HMDB_add_iso[ ,column_label] <= breaks_fwhm[length(breaks_fwhm)]), ]

  # remove adducts and isotopes, put internal standard at the beginning
  outlist_IS <- outlist[grep("IS", outlist[ , "CompoundName"], fixed=TRUE), ]
  outlist <- outlist[grep("HMDB", rownames(outlist), fixed=TRUE), ]
  outlist <- outlist[-grep("_", rownames(outlist), fixed=TRUE), ]
  outlist <- rbind(outlist_IS, outlist)
  # sort on m/z value
  outlist <- outlist[order(outlist[ ,column_label]), ]

  n <- dim(outlist)[1]
  # size of hmdb parts in lines:
  sub <- 2000
  end <- 0
  check <- 0

  # generate hmdb parts
  if (n >= sub & (floor(n/sub)) >= 2) {
    for (i in 1:floor(n/sub)){
      start <- -(sub-1) + i*sub
      end <- i*sub
      outlist_part <- outlist[c(start:end), ]
      save(outlist_part, file=paste0(scanmode, "_hmdb.", i, ".RData"))
    }
  }
}

# finish last hmdb part
start = end + 1
end = n

outlist_part = outlist[c(start:end),]
save(outlist_part, file=paste0(scanmode, "_hmdb.", i+1, ".RData"))
