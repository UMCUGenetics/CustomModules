#!/usr/bin/Rscript
# adapted from 6-peakGrouping.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

HMDB_part_file <- cmd_args[1]
SpecPeaks_file <- cmd_args[2]
pattern_file   <- cmd_args[3]
ppm <- as.numeric(cmd_args[4])

options(digits=16)

# load part of the HMDB
HMDB_add_iso <- get(load(HMDB_part_file)) 
# load(HMDB_part_file) 
# HMDB_add_iso <- outlist_part

# determine appropriate scanmode based on HMDB_part_file 
if (grepl("negative", basename(HMDB_part_file))) { scanmode <- "negative" } else
  if (grepl("positive", basename(HMDB_part_file))) { scanmode <- "positive" }

# determine batch number of HMDB part file
batch_number = strsplit(basename(HMDB_part_file), ".",fixed = TRUE)[[1]][2]

# load file with spectrum peaks
SpecPeaks_file <- paste0("SpectrumPeaks_", scanmode, ".RData")
load(SpecPeaks_file)
outlist.copy <- outlist.tot
rm(outlist.tot)

# load replication pattern
# load(paste0("./", scanmode, "_repl_pattern", ".RData"))
pattern_file <- paste0(scanmode, "_repl_pattern.RData")
load(pattern_file)
# load("./breaks.fwhm.RData")

# determine appropriate column name in HMDB part 
if (scanmode=="negative") { column_label <- "MNeg" } else { column_label <- "Mpos" }

# for debugging:
print(head(HMDB_add_iso))
print(scanmode)
print(column_label)

# Initialize
outpgrlist.identified <- NULL
outlist.grouped <- NULL

# First find peak groups identified based on HMDB masses
while (dim(HMDB_add_iso)[1] > 0) { 
  index <- 1
 
  # take one m/z value from the HMDB part and calculate mass tolerance
  reference_mass <- as.numeric(HMDB_add_iso[index, column_label])
  mass_tolerance <- (reference_mass * ppm) / 10^6
  
  # find the peaks in the dataset with corresponding m/z
  mzmed <- as.numeric(outlist.copy[ ,"mzmed.pkt"])
  selp <- which((mzmed > (reference_mass - mass_tolerance)) & (mzmed < (reference_mass + mass_tolerance)))
  tmplist <- outlist.copy[selp,,drop=FALSE]
  outlist.grouped <- rbind(outlist.grouped, tmplist)
  nrsamples <- length(selp)
  if (nrsamples > 0) {
    mzmed.pgrp <- mean(as.numeric(outlist.copy[selp, "mzmed.pkt"]))
    mzmin.pgrp <- reference_mass - mass_tolerance
    mzmax.pgrp <- reference_mass + mass_tolerance
    
    # determine fit quality fq
    fq.worst.pgrp <- as.numeric(max(outlist.copy[selp, "fq"]))
    fq.best.pgrp <- as.numeric(min(outlist.copy[selp, "fq"]))
    
    # set up object for intensities for all samples
    ints.allsamps <- rep(0, length(names(repl.pattern.filtered)))
    names(ints.allsamps) <- names(repl.pattern.filtered) 
    
    # Check for each sample if multiple peaks exist, if so take the sum of the intensities
    labels <- unique(tmplist[ ,"samplenr"])
    ints.allsamps[labels] <- as.vector(unlist(lapply(labels, function(x) { sum(as.numeric(tmplist[which(tmplist[ , "samplenr"]==x), "height.pkt"])) } )))
    
    # Initialize 
    assi_HMDB <- iso_HMDB <- HMDB_code <- NA
    tmplist.mass.iso <- tmplist.mass.adduct <- NULL
    
    # Identification: find all entries in HMDB part with mass within ppm range
    mass.all <- as.numeric(HMDB_add_iso[ , column_label])
    index <- which((mass.all > (reference_mass - mass_tolerance)) & (mass.all < (reference_mass + mass_tolerance)))
    tmplist.mass <- HMDB_add_iso[index,,drop=FALSE]
    
    if (dim(tmplist.mass)[1]>0) {
      # find isotope entries    
      index.iso <- grep(" iso ", tmplist.mass[, "CompoundName"], fixed = TRUE)
      if (length(index.iso) > 0){
        tmplist.mass.iso <- tmplist.mass[index.iso,,drop=FALSE]
        tmplist.mass <- tmplist.mass[-index.iso,,drop=FALSE]
      }
      
      if (dim(tmplist.mass)[1] > 0) {        
	# find adduct entries
        index.adduct <- grep(" [M", tmplist.mass[, "CompoundName"], fixed = TRUE)
        if (length(index.adduct) > 0) {
          tmplist.mass.adduct <- tmplist.mass[index.adduct,,drop=FALSE]
          tmplist.mass <- tmplist.mass[-index.adduct,,drop=FALSE]
        }  
      }  
      
      # Compose a list compounds, adducts or isotopes with corresponding m/z
      if (dim(tmplist.mass)[1]>0) {
        
        # metabolites
        assi_HMDB <- as.character(paste(as.character(tmplist.mass[, "CompoundName"]), collapse = ";"))
        HMDB_code <- as.character(paste(as.character(rownames(tmplist.mass)), collapse = ";"))
        theormz_HMDB <- as.numeric(tmplist.mass[1, column_label])
        
        # adducts of metabolites
        if (!is.null(tmplist.mass.adduct)) {
          if (dim(tmplist.mass.adduct)[1] > 0) {
            if (is.na(assi_HMDB)){
              assi_HMDB <- as.character(paste(as.character(tmplist.mass.adduct[, "CompoundName"]), collapse = ";"))
              HMDB_code <- as.character(paste(as.character(rownames(tmplist.mass.adduct)), collapse = ";"))
            } else {
              assi_HMDB <- paste(assi_HMDB, as.character(paste(as.character(tmplist.mass.adduct[, "CompoundName"]), collapse = ";")), sep = ";")
              HMDB_code <- paste(HMDB_code, as.character(paste(as.character(rownames(tmplist.mass.adduct)), collapse = ";")), sep = ";")
            }
	  }
	}
        
        # isotopes of metabolites
        if (!is.null(tmplist.mass.iso)) {
          if (dim(tmplist.mass.iso)[1]>0) {
            iso_HMDB <- as.character(paste(as.character(tmplist.mass.iso[, "CompoundName"]), collapse = ";"))
          }
	}
        
      # if no metabolites have the correct m/z, look for adducts and isotopes only
      } else if (!is.null(tmplist.mass.adduct)) {
        
        theormz_HMDB <- as.numeric(tmplist.mass.adduct[1, column_label])
        
        # adducts of metabolites
        if (!is.null(tmplist.mass.adduct)) {
          if (dim(tmplist.mass.adduct)[1] > 0) {
            if (is.na(assi_HMDB)) {
              assi_HMDB <- as.character(paste(as.character(tmplist.mass.adduct[ , "CompoundName"]), collapse = ";"))
              HMDB_code <- as.character(paste(as.character(rownames(tmplist.mass.adduct)), collapse = ";"))
            } else {
              assi_HMDB <- paste(assi_HMDB, as.character(paste(as.character(tmplist.mass.adduct[, "CompoundName"]), collapse = ";")), sep = ";")
              HMDB_code <- paste(HMDB_code, as.character(paste(as.character(rownames(tmplist.mass.adduct)), collapse = ";")), sep = ";")
            }
	  }
	}
        
        # isotopes of metabolites
        if (!is.null(tmplist.mass.iso)) {
          if (dim(tmplist.mass.iso)[1]>0) {
            iso_HMDB <- as.character(paste(as.character(tmplist.mass.iso[, "CompoundName"]), collapse = ";"))
          }
	}
        
      # if no metabolites or adducts can be found, only look for isotopes  
      } else if (!is.null(tmplist.mass.iso)) {
        
        if (dim(tmplist.mass.iso)[1]>0) {
          theormz_HMDB <- as.numeric(tmplist.mass.iso[1,column_label])
          iso_HMDB <- as.character(paste(as.character(tmplist.mass.iso[, "CompoundName"]), collapse = ";"))
        }
      }  
      
    }
   
    # combine all information 
    outpgrlist.identified <- rbind(outpgrlist.identified, cbind(data.frame(mzmed.pgrp, "fq.best"=fq.best.pgrp, "fq.worst"=fq.worst.pgrp, nrsamples, mzmin.pgrp, mzmax.pgrp),
                                         t(as.matrix(ints.allsamps)),
                                         data.frame(assi_HMDB, iso_HMDB, HMDB_code, theormz_HMDB)))
  }
 
  # remove index metabolite from HMDB part and continue while loop 
  HMDB_add_iso <- HMDB_add_iso[-index,]
  
}

# save peak list corresponding to masses in HMDB part
# save(outlist.grouped, file=paste0(batch_number, "_", scanmode, "_all.RData"))
# save peak group list, identified part
save(outpgrlist.identified, file=paste0(batch_number, "_", scanmode, "_identified.RData"))
