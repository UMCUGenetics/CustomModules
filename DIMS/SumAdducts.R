#!/usr/bin/Rscript
## adapted from 11-runSumAdducts.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

# collect_file <- cmd_args[1]
hmdbpart_main_file <- cmd_args[1]
scripts_dir <- cmd_args[2]
z_score <- as.numeric(cmd_args[3])
# outdir <- cmd_args[2]
# scanmode <- cmd_args[3]
# adducts <- cmd_args[4]

# for debugging:
print(hmdbpart_main_file)
# NB: scripts_dir not used yet, but function SumAdducts needs to be placed in AddOnFunctions folder
print(scripts_dir)
print(z_score)

if (grepl("positive_hmdb", hmdbpart_main_file)) { 
  scanmode <- "positive" 
  # for the adduct sum: include adducts M+Na (1) and M+K (2)
  adducts = c(1,2)
} else {
    if (grepl("negative_hmdb", hmdbpart_main_file)) { 
       scanmode <- "negative" 
       # for the adduct sum: include adduct M+Cl (1) 
       adducts <- c(1)
    }
}

# load input files
collect_file <- paste0("outlist_identified_", scanmode, ".RData")
load(collect_file)
repl_file <- paste0(scanmode, "_repl_pattern.RData")
load(repl_file)
outlist_part <- get(load(hmdbpart_main_file))

# adducts=as.vector(unlist(strsplit(adducts, ",",fixed = TRUE)))

# load(paste(outdir, "/outlist_identified_", scanmode, ".RData", sep=""))

# Local and on HPC
#batch = strsplit(file, "/",fixed = TRUE)[[1]]
#batch = batch[length(batch)]
#batch = strsplit(batch, ".",fixed = TRUE)[[1]][2]
batch_number = strsplit(basename(hmdbpart_main_file), ".",fixed = TRUE)[[1]][2]
print(batch_number)

outlist.tot <- unique(outlist.ident)

sumAdducts <- function(peaklist, theor_MZ, grpnames.long, adducts, batch_number, scanmode, outdir, z_score){
  #theor_MZ = outlist_part
  #grpnames.long = names(repl.pattern.filtered)
  #peaklist = outlist.ident
  #adducts = c(1) #for neg or c(1,2) for pos
  #batch <- 300
  #outdir <- "/Users/nunen/Documents/Metab/processed/zebrafish"
  #scanmode <- "negative"
  #z_score <- 0
  
  hmdb_codes <- rownames(theor_MZ)
  hmdb_names <- theor_MZ[,1, drop=FALSE]
  hmdb_names[] <- lapply(hmdb_names, as.character)
  
  # remove isotopes!!!
  index <- grep("HMDB", hmdb_codes, fixed=TRUE)
  hmdb_codes <- hmdb_codes[index]
  hmdb_names <- hmdb_names[index,]
  index = grep("_", rownames(hmdb_codes), fixed=TRUE)
  if (length(index) > 0) hmdb_codes <- hmdb_codes[-index]
  if (length(index) > 0) hmdb_names <- hmdb_names[-index]
  
  # negative
  names <- NULL
  adductsum <- NULL
  names_long <- NULL
  
  if (length(hmdb_codes) != 0) {
    
    for(i in 1:length(hmdb_codes)){
      #compound="HMDB00045"
      compound <- hmdb_codes[i]
      compound_plus <- c(compound,paste(compound, adducts, sep = "_"))
      
      # x=peaklist$HMDB_code[1]
      metab <- unlist(lapply(peaklist$HMDB_code, function(x) {(length(intersect(unlist(strsplit(as.vector(x),";")),compound_plus))>0)}))
      # peaklist[metab, "assi.hmdb"]
      # which(metab==TRUE)
      
      total <- c()
      
      # peaklist[metab, c("mzmed.pgrp", "HMDB_code", "C34.1")]
      # ints=peaklist[metab, c(7:(length(grpnames.long)+6))]
      if (z_score == 1) {
        ints <- peaklist[metab, c(15:(length(grpnames.long)+14))]
      } else {
        ints <- peaklist[metab, c(7:(length(grpnames.long)+6))]
      }
      total <- apply(ints, 2, sum)
      
      if (sum(total)!=0) {
        #message(i)
        names <- c(names, compound)
        adductsum <- rbind(adductsum,total)
        names_long <- c(names_long, hmdb_names[i])
      }
    }
    
    if (!is.null(adductsum)){ 
      rownames(adductsum) <- names
      adductsum <- cbind(adductsum, "HMDB_name"=names_long)
      save(adductsum, file = paste(scanmode, "_", batch_number, "_SummedAdducts.RData", sep=""))
    }
  }  
}


sumAdducts(outlist.tot, outlist_part, names(repl_pattern_filtered), adducts, batch_number, scanmode, outdir, z_score)

