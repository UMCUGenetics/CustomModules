#!/usr/bin/Rscript

# load required packages 
# none 

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n")

# outdir <- cmd_args[1]
# scanmode <- cmd_args[2]
scanmodes <- c("positive", "negative")

# Check if all jobs terminated correct!
notRun = NULL

# new: run for each scanmode
for (scanmode in scanmodes) {
  load(paste0("./repl.pattern.", scanmode, ".RData"))
  groupNames = names(repl.pattern.filtered)

  indir <- "./"
  object.files = list.files(indir, full.names=TRUE, pattern="*.RData")

  for (i in 1:length(groupNames)) {
    group <- paste0(indir, "/", paste0(paste(groupNames[i], scanmode, sep = "_"), ".RData"))
    if (!(group %in% object.files)) {
      notRun = c(notRun, group)
    }
  }

    #if (is.null(notRun)){
    cat("\nCollecting samples!")

    # negative
    filepath <- "./")
    files <- list.files(filepath, recursive=TRUE, full.names=TRUE, pattern=paste("*_", scanmode, ".RData",sep=""))

    outlist.tot=NULL
    for (i in 1:length(files)) {

      cat("\n", files[i])
      load(files[i])

      if (is.null(outlist.persample) || (dim(outlist.persample)[1]==0)){
        tmp=strsplit(files[i], "/")[[1]]
        fname = tmp[length(tmp)]
        #fname = strsplit(files[i], "/")[[1]][8]
        fname = strsplit(fname, ".RData")[[1]][1]
        fname = substr(fname, 13, nchar(fname))

        if (i == 1) { outlist.tot <- c(fname, rep("-1",5)) } else { outlist.tot <- rbind(outlist.tot, c(fname, rep("-1",5)))}
      } else {
        if (i == 1) { outlist.tot <- outlist.persample } else { outlist.tot <- rbind(outlist.tot, outlist.persample)}
      }
    }

    # remove negative values
    index=which(outlist.tot[,"height.pkt"]<=0)
    if (length(index)>0) outlist.tot = outlist.tot[-index,]
    index=which(outlist.tot[,"mzmed.pkt"]<=0)
    if (length(index)>0) outlist.tot = outlist.tot[-index,]

    # outdir_specpks <- paste(outdir, "5-specpks_all", sep = "/")
    # dir.create(outdir_specpks, showWarnings = F)
    # save(outlist.tot, file = paste(outdir_specpks, paste(scanmode, "RData", sep = "."), sep = "/"))
    save(outlist.tot, file = paste(scanmode, "RData", sep = ".")


    if (!is.null(notRun)){
      for (i in 1:length(notRun)){
        message(paste(notRun[i], "was not generated"))
      }
    }

}

