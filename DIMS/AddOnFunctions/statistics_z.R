statistics_z <- function(peaklist, sortCol, adducts){
  # peaklist=as.data.frame(outlist.adducts.HMDB)
  # plotdir="./results/plots/adducts"
  # filename="./results/allpgrps_stats.txt"
  # adducts=TRUE
  
  # peaklist=outlist.tot
  # sortCol=NULL
  # adducts=FALSE
  
  case_label = "P"
  control_label = "C"
  startcol = dim(peaklist)[2]+3 
  
  # calculate mean and sd for Control group
  ctrl.cols <- grep(control_label, colnames(peaklist),fixed = TRUE) # 5:41
  int.cols <- c(grep(control_label, colnames(peaklist),fixed = TRUE), grep(case_label, colnames(peaklist),fixed = TRUE))
  peaklist[,int.cols][peaklist[,int.cols]==0] = NA
  
  # tmp = data.matrix(peaklist[ , ctrl.cols], rownames.force = TRUE)
  tmp = peaklist[ , ctrl.cols, drop=FALSE]
  
  peaklist$avg.ctrls = apply(tmp, 1, function(x) mean(as.numeric(x),na.rm = TRUE))
  peaklist$sd.ctrls = apply(tmp, 1, function(x) sd(as.numeric(x),na.rm = TRUE))
  
  cnames.z = NULL
  
  for (i in int.cols) {
    # message(i)
    cname = colnames(peaklist)[i]
    cnames.z = c(cnames.z, paste(cname, "Zscore", sep="_"))
    zscores.1col = (as.numeric(as.vector(unlist(peaklist[ , i]))) - peaklist$avg.ctrls) / peaklist$sd.ctrls
    peaklist = cbind(peaklist, zscores.1col)
  }
  
  colnames(peaklist)[startcol:ncol(peaklist)] = cnames.z
  
  z.cols = grep("Zscore", colnames(peaklist),fixed = TRUE)
  
  if (!adducts){
    if ((dim(peaklist[, z.cols])[2]+6)!=(startcol-1)){
      ppmdev=array(1:dim(peaklist)[1], dim=c(dim(peaklist)[1]))
      
      # calculate ppm deviation
      for (i in 1:dim(peaklist)[1]){
        if (!is.na(peaklist$theormz_HMDB[i]) & !is.null(peaklist$theormz_HMDB[i]) & (peaklist$theormz_HMDB[i]!="")){
          ppmdev[i] = 10^6*(as.numeric(as.vector(peaklist$mzmed.pgrp[i]))-as.numeric(as.vector(peaklist$theormz_HMDB[i])))/as.numeric(as.vector(peaklist$theormz_HMDB[i]))
        } else {
          ppmdev[i]=NA
        }
      }
      
      peaklist = cbind(peaklist[, 1:6], ppmdev=ppmdev, peaklist[ , 7:ncol(peaklist)])
    }  
  }
  
  #peaklist = peaklist[order(peaklist[,sortCol]),]
  
  #   # Order on average Z-score
  #   tmp = peaklist[,grep("Zscore", colnames(peaklist))]
  #   tmp.p = tmp[,grep("P", colnames(tmp)),drop=FALSE]
  #   tmp.p.avg = apply(tmp.p, 1, mean)
  #   
  #   peaklist = cbind(peaklist, "avg.z.score"=tmp.p.avg)
  #   peaklist = peaklist[order(abs(tmp.p.avg), decreasing = TRUE),]
  
  return(peaklist)
  
}  
