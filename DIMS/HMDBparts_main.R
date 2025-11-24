# load package
load("argparse")

parser <- ArgumentParser(description = "HMDBparts_main")

parser$add_argument("--hmdb_db_file", dest = "hmdb_db_file",
                    help = "File path the HMDB databse file", required = TRUE)
parser$add_argument("--breaks_file", dest = "breaks_file",
                    help = "File containing all the breaks/bins", required = TRUE)

args <- parser$parse_args()

load(args$db_file)
load(args$breaks_file)

# Cut up HMDB minus adducts minus isotopes into small parts
scanmodes <- c("positive", "negative")
for (scanmode in scanmodes) {
  if (scanmode == "negative") {
    column_label <- "MNeg"
    hmdb_add_iso <- HMDB_add_iso.Neg
  } else if (scanmode == "positive") {
    column_label <- "Mpos"
    hmdb_add_iso <- HMDB_add_iso.Pos
  }

  # filter mass range measured
  outlist <- hmdb_add_iso[which(hmdb_add_iso[, column_label] >= breaks_fwhm[1] &
                                  hmdb_add_iso[, column_label] <= breaks_fwhm[length(breaks_fwhm)]), ]

  # remove adducts and isotopes, put internal standard at the beginning
  outlist <- outlist[grep("HMDB", rownames(outlist), fixed = TRUE), ]
  outlist <- outlist[-grep("_", rownames(outlist), fixed = TRUE), ]
  # sort on m/z value
  outlist <- outlist[order(outlist[, column_label]), ]
  nr_rows <- dim(outlist)[1]

  # size of hmdb parts in lines:
  sub <- 1000
  end <- 0
  check <- 0

  # generate hmdb parts
  if (nr_rows >= sub && (floor(nr_rows / sub)) >= 2) {
    for (i in 1:floor(nr_rows / sub)) {
      start <- -(sub - 1) + i * sub
      end <- i * sub
      outlist_part <- outlist[c(start:end), ]
      save(outlist_part, file = paste0(scanmode, "_hmdb_main.", i, ".RData"))
    }
  }

  # finish last hmdb part
  start <- end + 1
  end <- nr_rows

  outlist_part <- outlist[c(start:end), ]
  save(outlist_part, file = paste0(scanmode, "_hmdb_main.", i + 1, ".RData"))

}
