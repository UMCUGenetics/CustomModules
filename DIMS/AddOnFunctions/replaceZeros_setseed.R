replaceZeros <- function(file,scanmode,resol,outdir,thresh,scriptDir,ppm){
  # file="./results/grouping_rest/negative_1.RData"
  # scanmode= "negative"
  # scriptDir="./scripts"
  # resol=140000
  # thresh=2000
  # outdir="./results"
  
  control_label="C"
  
  source(paste(scriptDir, "AddOnFunctions/sourceDir.R", sep="/"))
  sourceDir(paste(scriptDir, "AddOnFunctions", sep="/"))
  
  dir.create(paste(outdir, "9-samplePeaksFilled", sep="/"), showWarnings = FALSE)
  
  # int.factor=1*10^5 # Number of x used to calc area under Gaussian (is not analytic)
  # scale=2 # Initial value used to estimate scaling parameter
  # width=1024
  # height=768
  
  # message(paste("file", file))
  # message(paste("scanmode", scanmode))
  # message(paste("resol", resol))
  # message(paste("outdir", outdir))
  # message(paste("thresh", thresh))
  # message(paste("scriptDir", scriptDir))
  
  load(paste0(outdir, "/repl.pattern.",scanmode, ".RData"))
  
  name = as.vector(unlist(strsplit(file, "/", fixed=TRUE)))
  name = name[length(name)]
  # message(paste("File name: ", name))
  
  # load samplePeaks
  # load  outpgrlist
  load(file)
  
  # #################################################################################
  # # filter on at least signal in two control samples
  # int.cols = grep(control_label, colnames(outpgrlist),fixed = TRUE)
  # # barplot(as.numeric(outpgrlist[753, int.cols]))
  # keep = NULL
  # keep = apply(outpgrlist, 1, function(x) if (length(which(as.numeric(x[int.cols]) > 0)) > 1) keep=c(keep,TRUE) else keep=c(keep,FALSE))
  # outpgrlist = outpgrlist[keep,]
  # #################################################################################
  
  ################################################################################
  # For now only replace zeros
  if (!is.null(outpgrlist)) {
    for (i in 1:length(names(repl.pattern.filtered))){
      samplePeaks=outpgrlist[,names(repl.pattern.filtered)[i]]
      index=which(samplePeaks<=0)
      if (!length(index)){
        next
      }
      for (j in 1:length(index)){
        area = generateGaussian(outpgrlist[index[j],"mzmed.pgrp"],thresh,resol,FALSE,scanmode,int.factor=1*10^5,1,1)$area
        # for testing purposes, add a fixed random seed
        set.seed(123)
        outpgrlist[index[j], names(repl.pattern.filtered)[i]] = rnorm(n=1, mean=area, sd=0.25*area)
      }
    }
  
  ################################################################################
  
  
  #################### identification #########################################################
    # load(paste(scriptDir, "../db/HMDB_add_iso_corrNaCl.RData", sep="/")) # E:\Metabolomics\LargeDataBase\Apr25_2016
    
    # Add average column
    outpgrlist = cbind(outpgrlist, "avg.int"=apply(outpgrlist[, 7:(ncol(outpgrlist)-4)], 1, mean))
    
    if (scanmode=="negative"){
      label = "MNeg"
      label2 = "Negative"
      # take out multiple NaCl adducts
      look4.add2 <- c("Cl", "Cl37", "For", "NaCl","KCl","H2PO4","HSO4","Na-H","K-H","H2O","I") # ,"min2H","min3H"
      # HMDB_add_iso=HMDB_add_iso.Neg
    } else {
      label = "Mpos"
      label2 = "Positive"
      # take out NaCl adducts
      look4.add2 <- c("Na", "K", "NaCl", "NH4","2Na-H","CH3OH","KCl","NaK-H") # ,"NaCl2","NaCl3","NaCl4","NaCl5")
      # HMDB_add_iso=HMDB_add_iso.Pos
    }
    
    # # Identification using large database
    # final.outlist.idpat = iden.code(outpgrlist, HMDB_add_iso, ppm=2, label)
    # message(paste(sum(final.outlist.idpat[ , "assi_HMDB"] != ""), "assigned peakgroups"))
    # message(paste(sum(final.outlist.idpat[ , "iso_HMDB"] != ""), "assigned isomeres"))
    
    # Identify noise peaks
    noise.MZ <- read.table(file="/hpc/dbg_mz/tools/db/TheoreticalMZ_NegPos_incNaCl.txt", sep="\t", header=TRUE, quote = "")
    noise.MZ <- noise.MZ[(noise.MZ[ , label] != 0), 1:4]
    
    # Replace "Negative" by "negative" in ident.hires.noise
    final.outlist.idpat2 = ident.hires.noise.HPC(outpgrlist, allAdducts, scanmode=label2, noise.MZ, look4=look4.add2, resol=resol, slope=0, incpt=0, ppm.fixed=ppm, ppm.iso.fixed=ppm)
    # message(paste(sum(final.outlist.idpat2[ , "assi"] != ""), "assigned noise peaks"))
    tmp <- final.outlist.idpat2[ , c("assi", "theormz")]
    colnames(tmp) <- c("assi_noise",  "theormz_noise")
    
    final.outlist.idpat3 <- cbind(outpgrlist, tmp)
    #############################################################################################
    
    # message(paste("File saved: ", paste(outdir, "/samplePeaksFilled/", name, sep="")))
    save(final.outlist.idpat3, file=paste(outdir, "/9-samplePeaksFilled/", name, sep=""))
  }
}
