#!/usr/bin/Rscript
## adapted from 8-peakGrouping.rest.R

# load required packages 
# none 

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)

resol <- as.numeric(cmd_args[1])
ppm <- as.numeric(cmd_args[2])
outdir <- "./"
options(digits=16)

print(list.files(outdir, pattern="RData"))

scanmodes <- c("positive", "negative")

# function for grouping unidentified peaks
groupingRest <- function(outdir, unidentified_peaklist, scanmode, ppm) {
  outlist.copy <- get(load(unidentified_peaklist))
  # batch = strsplit(unidentified_peaklist, ".",fixed = TRUE)[[1]][2] 
  load(paste0("./", scanmode, "_repl_pattern.RData"))
  
  outpgrlist <- NULL
  
  # Then group on highest peaks
  range <- ppm*1e-06
  startcol <- 7
  
  # while (max(as.numeric(outlist.copy[ , "height.pkt"])) > 0 ) {
  while (dim(outlist.copy)[1] > 0) {
    
    sel <- which(as.numeric(outlist.copy[ , "height.pkt"]) == max(as.numeric(outlist.copy[ , "height.pkt"])))[1]
    
    # ppm range around max
    mzref <- as.numeric(outlist.copy[sel, "mzmed.pkt"])
    pkmin <- -(range*mzref - mzref)
    pkmax <- 2*mzref-pkmin
    
    selp <- as.numeric(outlist.copy[ , "mzmed.pkt"]) > pkmin & as.numeric(outlist.copy[ , "mzmed.pkt"]) < pkmax
    tmplist <- outlist.copy[selp,,drop=FALSE]
    
    nrsamples <- length(unique(tmplist[,"samplenr"]))
    if (nrsamples > 0) {
      
      mzmed.pgrp <- mean(as.numeric(outlist.copy[selp, "mzmed.pkt"]))
      mzmin.pgrp <- -(range*mzmed.pgrp - mzmed.pgrp)
      mzmax.pgrp <- 2*mzmed.pgrp - mzmin.pgrp
      
      selp <- as.numeric(outlist.copy[ , "mzmed.pkt"]) > mzmin.pgrp & as.numeric(outlist.copy[ , "mzmed.pkt"]) < mzmax.pgrp
      tmplist <- outlist.copy[selp,,drop=FALSE]
      
      # remove used peaks!!!
      tmp <- as.vector(which(tmplist[,"height.pkt"]==-1))
      if (length(tmp)>0) tmplist<-tmplist[-tmp,,drop=FALSE]
      
      nrsamples <- length(unique(tmplist[,"samplenr"]))
      
      fq.worst.pgrp <- as.numeric(max(outlist.copy[selp, "fq"]))
      fq.best.pgrp <- as.numeric(min(outlist.copy[selp, "fq"]))
      ints.allsamps <- rep(0, length(names(repl_pattern_filtered)))
      names(ints.allsamps) <- names(repl_pattern_filtered) # same order as sample list!!!
      
      # Check for each sample if multiple peaks exists, if so take the sum!
      labels <- unique(tmplist[,"samplenr"])
      ints.allsamps[labels] <- as.vector(unlist(lapply(labels, function(x) {sum(as.numeric(tmplist[which(tmplist[ , "samplenr"]==x), "height.pkt"]))})))
      
      outpgrlist <- rbind(outpgrlist, c(mzmed.pgrp, fq.best.pgrp, fq.worst.pgrp, nrsamples, mzmin.pgrp, mzmax.pgrp, ints.allsamps,NA,NA,NA,NA))
    }
   
    outlist.copy <- outlist.copy[-which(selp==TRUE),,drop=FALSE]
  }
  
  outpgrlist <- as.data.frame(outpgrlist)  # ignore warnings of duplicate row names
  colnames(outpgrlist)[1:6] <- c("mzmed.pgrp", "fq.best", "fq.worst", "nrsamples", "mzmin.pgrp", "mzmax.pgrp")
  colnames(outpgrlist)[(length(repl_pattern_filtered)+7):ncol(outpgrlist)] <- c("assi_HMDB", "iso_HMDB", "HMDB_code", "theormz_HMDB")
 
  return(outpgrlist) 
  
}


for (scanmode in scanmodes) {
  unidentified_peaklist <- paste0("./SpectrumPeaks_", scanmode, "_Unidentified.RData") 
  # generate peak group lists of the unidentified peaks
  outpgrlist <- groupingRest(outdir, unidentified_peaklist, scanmode, ppm=ppm)

  # save output in RData format for further processing  
  save(outpgrlist, file=paste0("PeakGroupList_", scanmode, "_Unidentified.RData"))
  write.table(outpgrlist, file=paste0("PeakGroupList_", scanmode, "_Unidentified.txt"))
}

