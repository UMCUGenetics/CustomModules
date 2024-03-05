generateExcelFile <- function(peaklist,  fileName, sub, plot = TRUE) {
  # plotdir=file.path(plotdir)
  # imageNum=2
  # fileName="./results/xls/Pos_allpgrps_identified"
  # subName=c("","_box")
  
  end <- 0
  i <- 0
  
  if (dim(peaklist)[1]>=sub & (sub>0)){
    for (i in 1:floor(dim(peaklist)[1]/sub)){
      start=-(sub-1)+i*sub
      end=i*sub
      message(paste0(start, ":", end))
      
      genExcelFileV3(peaklist[c(start:end),], paste(fileName, sep="_"), plot)
    }
  }  
  start = end + 1
  end = dim(peaklist)[1]
  message(start)
  message(end)
  genExcelFileV3(peaklist[c(start:end),], paste(fileName, i+1, sep="_"), plot)
}
