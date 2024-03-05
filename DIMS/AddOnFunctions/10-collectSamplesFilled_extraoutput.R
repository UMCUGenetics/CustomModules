#!/usr/bin/Rscript

.libPaths(new = "/hpc/local/CentOS7/dbg_mz/R_libs/3.2.2")

# load required packages 
# none 

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n")

outdir <- cmd_args[1]
scanmode <- cmd_args[2]
normalization <- cmd_args[3]
scripts <- cmd_args[4]
z_score <- as.numeric(cmd_args[5])
ppm <- as.numeric(cmd_args[6])

#outdir <- "/Users/nunen/Documents/Metab/processed/test_old"
#scanmode <- "negative"
#normalization <- "disabled"
#scripts <- "/Users/nunen/Documents/Metab/DIMS/scripts"
#db <- "/Users/nunen/Documents/Metab/DIMS/db/HMDB_add_iso_corrNaCl_withIS_withC5OH.RData"
#z_score <- 0

object.files = list.files(paste(outdir, "9-samplePeaksFilled", sep="/"), full.names=TRUE, pattern=scanmode)
outlist.tot=NULL
for (i in 1:length(object.files)) {
  load(object.files[i])
  print(print(object.files[i]))
  outlist.tot = rbind(outlist.tot, final.outlist.idpat3)
}

source(paste(scripts, "AddOnFunctions/sourceDir.R", sep="/"))
sourceDir(paste(scripts, "AddOnFunctions", sep="/"))

# remove duplicates
outlist.tot = mergeDuplicatedRows(outlist.tot)

# sort on mass
outlist.tot = outlist.tot[order(outlist.tot[,"mzmed.pgrp"]),]

# normalization
load(paste0(outdir, "/repl.pattern.",scanmode,".RData"))

if (normalization != "disabled") {
  outlist.tot = normalization_2.1(outlist.tot, fileName, names(repl.pattern.filtered), on=normalization, assi_label="assi_HMDB")
}

if (z_score == 1) {
  outlist.stats = statistics_z(outlist.tot, sortCol=NULL, adducts=FALSE)
  nr.removed.samples=length(which(repl.pattern.filtered[]=="character(0)"))
  order.index.int=order(colnames(outlist.stats)[8:(length(repl.pattern.filtered)-nr.removed.samples+7)])
  outlist.stats.more = cbind(outlist.stats[,1:7],
                             outlist.stats[,(length(repl.pattern.filtered)-nr.removed.samples+8):(length(repl.pattern.filtered)-nr.removed.samples+8+6)],
                             outlist.stats[,8:(length(repl.pattern.filtered)-nr.removed.samples+7)][order.index.int],
                             outlist.stats[,(length(repl.pattern.filtered)-nr.removed.samples+5+10):ncol(outlist.stats)])
  
  tmp.index=grep("_Zscore", colnames(outlist.stats.more), fixed = TRUE)
  tmp.index.order=order(colnames(outlist.stats.more[,tmp.index]))
  tmp = outlist.stats.more[,tmp.index[tmp.index.order]]
  outlist.stats.more=outlist.stats.more[,-tmp.index]
  outlist.stats.more=cbind(outlist.stats.more,tmp)
  outlist.tot = outlist.stats.more
}

# filter identified compounds
index.1=which((outlist.tot[,"assi_HMDB"]!="") & (!is.na(outlist.tot[,"assi_HMDB"])))
index.2=which((outlist.tot[,"iso_HMDB"]!="") & (!is.na(outlist.tot[,"iso_HMDB"])))
index=union(index.1,index.2)
outlist.ident = outlist.tot[index,]
outlist.not.ident = outlist.tot[-index,]

if (z_score == 1) {
  outlist.ident$ppmdev=as.numeric(outlist.ident$ppmdev)
  outlist.ident <- outlist.ident[which(outlist.ident["ppmdev"] >= -ppm & outlist.ident["ppmdev"] <= ppm),]
}
# NAs in theormz_noise <======================================================================= uitzoeken!!!
outlist.ident$theormz_noise[which(is.na(outlist.ident$theormz_noise))] = 0
outlist.ident$theormz_noise=as.numeric(outlist.ident$theormz_noise)
outlist.ident$theormz_noise[which(is.na(outlist.ident$theormz_noise))] = 0
outlist.ident$theormz_noise=as.numeric(outlist.ident$theormz_noise)

save(outlist.not.ident, outlist.ident, file=paste(outdir, "/outlist_identified_", scanmode, ".RData", sep=""))

# Extra output in Excel-readable format:
remove_columns <- c("fq.best", "fq.worst", "mzmin.pgrp", "mzmax.pgrp")
remove_colindex <- which(colnames(outlist.ident) %in% remove_columns)
outlist.ident <- outlist.ident[ , -remove_colindex]
write.table(outlist.ident, file=paste0(outdir, "/outlist_identified_", scanmode, ".txt"), sep="\t", row.names = FALSE)
remove_colindex <- which(colnames(outlist.not.ident) %in% remove_columns)
outlist.not.ident <- outlist.not.ident[ , -remove_colindex]
write.table(outlist.not.ident, file=paste0(outdir, "/outlist_not_identified_", scanmode, ".txt"), sep="\t", row.names = FALSE)

