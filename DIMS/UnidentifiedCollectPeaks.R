#!/usr/bin/Rscript
## adapted from 7-collectSamplesGroupedHMDB.R

# load required packages 
# none 

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)

ppm <- as.numeric(cmd_args[1])
outdir <- "./" 

scanmodes <- c("positive", "negative")

for (scanmode in scanmodes) {
  # get list of all files that contain lists of peaks that were used in identified peak grouping
  files <- list.files("./", pattern = paste(scanmode, "_peaks_used.RData", sep=""))
  # load the list of all peaks
  load(paste0("SpectrumPeaks_", scanmode, ".RData")) #outlist.tot

  # Make a list of indexes of peaks that have been identified, then remove these from the peaklist.
  remove <- NULL
  for (i in 1:length(files)) {
    load(files[i]) # outlist.grouped, now called list_of_peaks_used_in_peak_groups_identified
    remove <- c(remove, which(outlist.tot[ ,"mzmed.pkt"] %in% list_of_peaks_used_in_peak_groups_identified[i ,"mzmed.pkt"]))
  }
  outlist.rest <- outlist.tot[-remove, ]

  # sort on mass
  outlist <- outlist.rest[order(as.numeric(outlist.rest[ ,"mzmed.pkt"])),]

  save(outlist, file=paste0(outdir, "/SpectrumPeaks_", scanmode, "_Unidentified.RData"))

  # cut the unidentified peak list in parts of part_len length for parallel processing
  # NB: the while statement below gives an error message. Skip the cutting into parts.
  # part_len <- 10000
  # nr_peaks <- dim(outlist)[1]
  # end <- 0
  # min_1_last <- part_len
  # check <- 0
  # outlist_i_min_1 <- NULL
  # part_nr <- 0

  # if (nr_peaks >= part_len & (floor(nr_peaks/part_len) - 1) >= 2) {
  #   for (part_nr in 2:floor(nr_peaks/part_len) -1 ) {
  #     start <- -(part_len-1) + part_nr*part_len
  #     end <- part_nr*part_len
  #     
  #     if (part_nr > 1) {
  #       outlist_i <- outlist[c(start:end),]
  # 
  #         n_moved <- 0
       
  #      # Calculate ppm and replace border, avoid cut within peakgroup!
  #      while ((as.numeric(outlist_i[1, "mzmed.pkt"]) - as.numeric(outlist_i_min_1[min_1_last, "mzmed.pkt"])) * 1e+06/as.numeric(outlist_i[1, "mzmed.pkt"]) < ppm) {
  #        outlist_i_min_1 <- rbind(outlist_i_min_1, outlist_i[1,])
  #        outlist_i <- outlist_i[-1, ]
  #        n_moved <- n_moved + 1
  #      }
      
  #      # message(paste("Process", i-1,":", dim(outlist_i_min_1)[1]))
  #      save(outlist_i_min_1, file = paste(outdir, paste(scanmode, paste("outlist_i_min_1", part_nr-1, "RData", sep="."), sep="_"), sep="/"))
  #      check <- check + dim(outlist_i_min_1)[1]
      
  #      outlist_i_min_1 <- outlist_i
  #      min_1_last <- dim(outlist_i_min_1)[1]
      
  #    } else {
  #      outlist_i_min_1 = outlist[c(start:end), ]
  #    }
  #  }
  #}

  #start <- end + 1
  #end <- nr_peaks
  #outlist_i <- outlist[c(start:end), ]
  #n_moved <- 0

  #if(!is.null(outlist_i_min_1)){
  #  # Calculate ppm and replace border, avoid cut within peakgroup!
  #  while ((as.numeric(outlist_i[1,"mzmed.pkt"]) - as.numeric(outlist_i_min_1[min_1_last,"mzmed.pkt"]))*1e+06/as.numeric(outlist_i[1,"mzmed.pkt"]) < ppm) {
  #    outlist_i_min_1 = rbind(outlist_i_min_1, outlist_i[1,])
  #    outlist_i = outlist_i[-1,]
  #    n_moved = n_moved + 1
  #  }
  
  #  cat(paste("Process", i+1-1,":", dim(outlist_i_min_1)[1]))
  #  save(outlist_i_min_1, file=paste(outdir, paste(scanmode, paste("outlist_i_min_1",i,"RData", sep="."), sep="_"), sep="/"))
  #  check=check+dim(outlist_i_min_1)[1]
  #}

  #outlist_i_min_1=outlist_i
  #cat("Process", i+2-1,":", dim(outlist_i_min_1)[1], "\n")
  #save(outlist_i_min_1, file=paste(outdir, paste(scanmode, paste("outlist_i_min_1",i+1,"RData", sep="."), sep="_"), sep="/"))
  #}

  #check <- check + dim(outlist_i_min_1)[1]
  #if (check==dim(outlist)[1]){
  #  cat("Check is oke!\n")
  #} else {
  #  cat("Check is failed!\n")
}
