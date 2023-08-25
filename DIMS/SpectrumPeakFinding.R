#!/usr/bin/Rscript
# adapted from 5-collectSamples.R 

# define parameters 
scanmodes <- c("positive", "negative")

# Check whether all jobs terminated correct!
notRun = NULL

print("one")

# collect spectrum peaks for each scanmode
for (scanmode in scanmodes) {
  print(scanmode)
  # load peak lists of all biological samples
  input_dir <- getwd() # "./"
  peaklist_files = list.files(input_dir, full.names=TRUE, pattern=paste("*_", scanmode, ".RData",sep=""))

  # get sample names
  load(paste0("./", scanmode, "_repl_pattern", ".RData"))
  groupNames = names(repl_pattern_filtered)
  for (i in 1:length(groupNames)) {
    group <- paste0(input_dir, "/", paste0(paste(groupNames[i], scanmode, sep = "_"), ".RData"))
    if (!(group %in% peaklist_files)) {
      notRun = c(notRun, group)
    }
  }
  
  cat("\nCollecting samples")

  outlist.tot=NULL
  for (i in 1:length(peaklist_files)) {
    cat("\n", peaklist_files[i])
    load(peaklist_files[i])
    if (is.null(outlist.persample) || (dim(outlist.persample)[1] == 0)) {
      tmp=strsplit(peaklist_files[i], "/")[[1]]
      fname = tmp[length(tmp)]
      fname = strsplit(fname, ".RData")[[1]][1]
      fname = substr(fname, 13, nchar(fname))
      if (i == 1) { 
        outlist.tot <- c(fname, rep("-1",5)) 
      } else { 
	outlist.tot <- rbind(outlist.tot, c(fname, rep("-1",5)))
      }
    } else {
      if (i == 1) { 
	outlist.tot <- outlist.persample 
      } else { 
	outlist.tot <- rbind(outlist.tot, outlist.persample) 
      }
    }
  }
  
  # remove negative values
  index <- which(outlist.tot[ ,"height.pkt"] <= 0)
  if (length(index) > 0) outlist.tot <- outlist.tot[-index, ]
  index <- which(outlist.tot[ ,"mzmed.pkt"] <= 0)
  if (length(index) > 0) outlist.tot <- outlist.tot[-index, ]

  save(outlist.tot, file = paste0("./SpectrumPeaks_", scanmode, ".RData"))

  if (!is.null(notRun)){
    for (i in 1:length(notRun)){
      message(paste(notRun[i], "was not generated"))
    }
  }
}

