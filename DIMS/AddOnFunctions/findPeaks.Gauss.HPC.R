### fit Gaussian estimate mean and integrate to obtain intensity
findPeaks.Gauss.HPC <- function(plist, breaks.fwhm, int.factor, scale, resol, outdir, scanmode, plot, thresh, width, height) {
  sampname <- colnames(plist)[1]
  
  range <- as.vector(plist)
  names(range) <- rownames(plist)
  
  values <- list("mean"=NULL, "area"=NULL, "nr"=NULL, "min"=NULL, "max"=NULL, "qual"=NULL, "spikes"=0)
  
  values <- searchMZRange(range, values, int.factor, scale, resol, outdir, sampname, scanmode, plot, width, height, thresh)
  
  outlist.persample <- NULL
  outlist.persample <- cbind("samplenr"=values$nr, "mzmed.pkt"=values$mean, "fq"=values$qual, "mzmin.pkt"=values$min, "mzmax.pkt"=values$max, "height.pkt"=values$area)

  index <- which(outlist.persample[ ,"height.pkt"]==0)
  if (length(index) > 0) {
    outlist.persample <- outlist.persample[-index,]
  }
  
  # save(outlist.persample, file=paste(outdir, paste(sampname, "_", scanmode, ".RData", sep=""), sep="/"))
  save(outlist.persample, file=paste("./", sampname, "_", scanmode, ".RData", sep=""))
  
  cat(paste("There were", values$spikes, "spikes!"))
}
