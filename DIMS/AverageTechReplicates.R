#!/usr/bin/Rscript
# adapted from 3-AverageTechReplicates.R

# load required packages 
# none 

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)
for (arg in cmd_args) cat("  ", arg, "\n", sep="")

init_filepath <- cmd_args[1]
nr_replicates <- as.numeric(cmd_args[2])
thresh2remove <- 2000
dimsThresh <- 100

removeFromRepl.pat <- function(bad_samples, repl.pattern, nr_replicates) {
  # bad_samples=remove_pos
  
  tmp = repl.pattern
  
  removeFromGroup=NULL
  
  for (i in 1:length(tmp)){
    tmp2 = repl.pattern[[i]]
    
    remove=NULL
    
    for (j in 1:length(tmp2)){
      if (tmp2[j] %in% bad_samples){
        #cat(tmp2[j])
        #cat(paste("remove",tmp2[j]))
        #cat(paste("remove i",i))
        #cat(paste("remove j",j))
        remove = c(remove, j)
      }
    }
    
    if (length(remove)==nr_replicates) removeFromGroup=c(removeFromGroup,i)
    if (!is.null(remove)) repl.pattern[[i]]=repl.pattern[[i]][-remove]
  }
  
  if (length(removeFromGroup)!=0) {
    repl.pattern=repl.pattern[-removeFromGroup]
  }
  
  return(list("pattern"=repl.pattern))
}


# get repl.pattern
load("./init.RData")

remove_neg=NULL
remove_pos=NULL
cat("Pklist sum threshold to remove technical replicate:", thresh2remove, "\n")
for (i in 1:length(repl.pattern)) {
  techRepsArray.pos = NULL
  techRepsArray.neg = NULL
  
  tech_reps = as.vector(unlist(repl.pattern[i]))
  sum_neg=0
  sum_pos=0
  n_pos=0
  n_neg=0
  cat("\n\nNow sample ", i, " from replication pattern with length ", length(repl.pattern))
  for (j in 1:length(tech_reps)) {
    load(paste("./", tech_reps[j], ".RData", sep=""))
    cat("\n\nParsing", tech_reps[j])

    cat("\n\tNegative pklist sum",sum(pklist$neg[,1]))
    if (sum(pklist$neg[,1])<thresh2remove){
      cat(" ... Removed")
      remove_neg=c(remove_neg, tech_reps[j])
    } else {
      n_neg=n_neg+1
      sum_neg=sum_neg+pklist$neg
    }
    
    techRepsArray.neg = cbind(techRepsArray.neg, pklist$neg)
    
    cat("\n\tPositive pklist sum",sum(pklist$pos[,1]))
    if (sum(pklist$pos[,1])<thresh2remove){
      cat(" ... Removed")
      remove_pos=c(remove_pos, tech_reps[j])
    } else {
      n_pos=n_pos+1
      sum_pos=sum_pos+pklist$pos
    }
    
    techRepsArray.pos = cbind(techRepsArray.pos, pklist$pos)
  }
  
  # filter within bins on at least signal in more than one tech. rep.!!!
  if (!is.null(dim(sum_pos))) sum_pos[apply(techRepsArray.pos,1,function(x) length(which(x>dimsThresh))==1),1]=0
  if (!is.null(dim(sum_neg))) sum_neg[apply(techRepsArray.neg,1,function(x) length(which(x>dimsThresh))==1),1]=0
  
  if (n_neg != 0){
    sum_neg[,1] <- sum_neg[,1]/n_neg
    colnames(sum_neg) <- names(repl.pattern)[i]
    save(sum_neg, file=paste("./", names(repl.pattern)[i], "_neg.RData", sep=""))
  }
  if (n_pos != 0) {
    sum_pos[,1] <- sum_pos[,1]/n_pos
    colnames(sum_pos) <- names(repl.pattern)[i]
    save(sum_pos, file=paste("./", names(repl.pattern)[i], "_pos.RData", sep=""))
  }
}

retVal <- removeFromRepl.pat(remove_neg, repl.pattern, nr_replicates)
repl.pattern.filtered <- retVal$pattern
save(repl.pattern.filtered, file="./repl_pattern_negative.RData")
write.table(remove_neg, file="./miss_infusions_neg.txt", row.names=FALSE, col.names=FALSE ,sep= "\t")

retVal <- removeFromRepl.pat(remove_pos, repl.pattern, nr_replicates)
repl.pattern.filtered <- retVal$pattern
save(repl.pattern.filtered, file="./repl_pattern_positive.RData")
write.table(remove_pos, file="./miss_infusions_pos.txt", row.names=FALSE, col.names=FALSE ,sep= "\t")


