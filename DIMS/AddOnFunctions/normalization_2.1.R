# normalization_2.1(outlist.pos.id, "Intensity_all_peaks_positive_norm", groupNames.pos, on="total", assi_label="assi_HMDB")

normalization_2.1 <- function(data, filename, groupNames, on="total", assi_label="assi"){
  # data=outlist.pos.id
  # filename = "Intensity_all_peaks_positive_norm"
  # groupNames = groupNames.pos
  # on="total"
  # assi_label="assi_HMDB"
  
  lastcol = length(groupNames) + 6
  before = data[ ,c(7:lastcol)]
  
  if (on=="total_IS") {
    
    assi = which(colnames(data)==assi_label)
    data.int <- data[ ,c(assi,7:lastcol)]  # assi and samples columns
    data.assi = data.int[grep("(IS", data.int[,1], ignore.case=FALSE, fixed = TRUE),]
    
  } else if (on=="total_ident"){
    
    assi = which(colnames(data)==assi_label)
    assi.hmdb = which(colnames(data)=="assi.hmdb")
    index = sort(union(which(data[,assi]!=""), which(data[,assi.hmdb]!="")))
    
    data.int = data[ ,c(assi,7:lastcol)]  # assi and samples columns
    data.assi = data.int[index,]
    
  } else if (on=="total") {
    
    assi = which(colnames(data)==assi_label)
    data.int = data[ ,c(assi,7:lastcol)]  # assi and samples columns
    data.assi = data.int
    
  }
  
  sum <- 0
  for (c in 2:ncol(data.assi)) {
    sum <- sum + sum(as.numeric(data.assi[,c]))
  }
  average <- sum/(ncol(data.assi)-1)
  for (c in 2:ncol(data.int)) {
    factor <- sum(as.numeric(data.assi[,c]))/average
    if (factor==0) {
      data.int[ ,c]=0
      cat(colnames(data.int)[c])
      cat("factor==0 !!!")
    } else {
      data.int[ ,c] <- as.numeric(data.int[ ,c])/factor
    }
  }
  
  # colnames(data.int[,2:ncol(data.int)])
  
  if (dim(data)[2]==lastcol){
    final.outlist.Pos.idpat.norm <- cbind(data[,1:6],data.int[,2:ncol(data.int)])
  } else {
    final.outlist.Pos.idpat.norm <- cbind(data[,1:6],data.int[,2:ncol(data.int)],data[,(lastcol + 1):ncol(data)])
  }
  
  #outdir="./results/normalization"
  #dir.create(outdir, showWarnings = FALSE)
  
  #CairoPNG(filename=paste(outdir, paste(filename, "_before.png", sep=""), sep="/"), width, height)
  #sub=apply(before,2, function(x) sum(as.numeric(x)))
  #barplot(as.vector(unlist(sub)), main="Not normalized",names.arg = colnames(before),las=2)
  #dev.off()
  
  #CairoPNG(filename=paste(outdir, paste(filename, "_", on, ".png", sep=""), sep="/"), width, height)
  #sub=apply(final.outlist.Pos.idpat.norm[,c(7:lastcol)],2, function(x) sum(as.numeric(x)))
  #barplot(as.vector(unlist(sub)), main=filename,names.arg = colnames(final.outlist.Pos.idpat.norm[,c(7:lastcol)]),las=2)
  #dev.off()
  
  return(final.outlist.Pos.idpat.norm)
  
}
