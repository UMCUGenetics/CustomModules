#!/usr/bin/Rscript
# adapted from 3-AverageTechReplicates.R

# define parameters 
cmd_args <- commandArgs(trailingOnly = TRUE)

init_filepath <- cmd_args[1]
nr_replicates <- as.numeric(cmd_args[2])
thresh2remove <- 1000000000
dimsThresh    <- 100
run_name      <- cmd_args[3]
dims_matrix   <- cmd_args[4]
highest_mz    <- as.numeric(cmd_args[5])

print(run_name)
print(dims_matrix)
print(highest_mz)

if (dims_matrix == "DBS") {
    thresh2remove <- 50000000
}
if (highest_mz > 700) {
    thresh2remove <- 1000000
}

library("ggplot2")
library("gridExtra")

removeFromRepl.pat <- function(bad_samples, repl_pattern, nr_replicates) {
  # bad_samples=remove_pos
  tmp <- repl_pattern
  removeFromGroup <- NULL
  
  for (i in 1:length(tmp)){
    tmp2 <- repl_pattern[[i]]
    remove <- NULL
    for (j in 1:length(tmp2)){
      if (tmp2[j] %in% bad_samples){
        remove = c(remove, j)
      }
    }
    
    if (length(remove) == nr_replicates) removeFromGroup <- c(removeFromGroup,i)
    if (!is.null(remove)) repl_pattern[[i]] <- repl_pattern[[i]][-remove]
  }
  
  if (length(removeFromGroup)!=0) {
    repl_pattern <- repl_pattern[-removeFromGroup]
  }
  
  return(list("pattern"=repl_pattern))
}


# get repl_pattern
load("./init.RData")

remove_neg <- NULL
remove_pos <- NULL
cat("Pklist sum threshold to remove technical replicate:", thresh2remove, "\n")
for (i in 1:length(repl_pattern)) {
  techRepsArray.pos <- NULL
  techRepsArray.neg <- NULL
  
  tech_reps <- as.vector(unlist(repl_pattern[i]))
  sum_neg <- 0
  sum_pos <- 0
  n_pos <- 0
  n_neg <- 0
  for (j in 1:length(tech_reps)) {
    load(paste("./", tech_reps[j], ".RData", sep=""))
    cat("\n\tNegative peak_list sum", sum(peak_list$neg[,1]))
    if (sum(peak_list$neg[,1]) < thresh2remove){
      cat(" ... Removed")
      remove_neg <- c(remove_neg, tech_reps[j])
    } else {
      n_neg <- n_neg + 1
      sum_neg <- sum_neg + peak_list$neg
    }
    
    techRepsArray.neg <- cbind(techRepsArray.neg, peak_list$neg)
    
    cat("\n\tPositive peak_list sum", sum(peak_list$pos[,1]))
    if (sum(peak_list$pos[,1]) < thresh2remove){
      cat(" ... Removed")
      remove_pos <- c(remove_pos, tech_reps[j])
    } else {
      n_pos <- n_pos + 1
      sum_pos <- sum_pos + peak_list$pos
    }
    
    techRepsArray.pos <- cbind(techRepsArray.pos, peak_list$pos)
  }
  
  # filter within bins on at least signal in more than one tech. rep.!!!
  if (!is.null(dim(sum_pos))) sum_pos[apply(techRepsArray.pos,1,function(x) length(which(x>dimsThresh))==1),1]=0
  if (!is.null(dim(sum_neg))) sum_neg[apply(techRepsArray.neg,1,function(x) length(which(x>dimsThresh))==1),1]=0
  
  if (n_neg != 0){
    sum_neg[,1] <- sum_neg[,1]/n_neg
    colnames(sum_neg) <- names(repl_pattern)[i]
    save(sum_neg, file=paste("./", names(repl_pattern)[i], "_neg_avg.RData", sep=""))
  }
  if (n_pos != 0) {
    sum_pos[,1] <- sum_pos[,1]/n_pos
    colnames(sum_pos) <- names(repl_pattern)[i]
    save(sum_pos, file=paste("./", names(repl_pattern)[i], "_pos_avg.RData", sep=""))
  }
}

retVal <- removeFromRepl.pat(remove_neg, repl_pattern, nr_replicates)
repl_pattern_filtered <- retVal$pattern
save(repl_pattern_filtered, file = "./negative_repl_pattern.RData")
write.table(remove_neg, file = "./miss_infusions_negative.txt", row.names=FALSE, col.names=FALSE , sep= "\t")

retVal <- removeFromRepl.pat(remove_pos, repl_pattern, nr_replicates)
repl_pattern_filtered <- retVal$pattern
save(repl_pattern_filtered, file = "./positive_repl_pattern.RData")
write.table(remove_pos, file = "./miss_infusions_positive.txt", row.names=FALSE, col.names=FALSE , sep= "\t")

# New: generate TIC plots
# get all txt files
TIC_files = list.files("./", full.names=TRUE, pattern="*TIC.txt")
all_samps <- sub('_TIC\\..*$', '', basename(TIC_files))

# determine maximum intensity
highest_tic_max <- 0
for (file in TIC_files) {
  tic <- read.table(file)
  this_tic_max <- max(tic$TIC)
  if (this_tic_max > highest_tic_max) {
    highest_tic_max <- this_tic_max
    max_sample <- sub('_TIC\\..*$', '', basename(file))
    }
}

tic_plot_list <- list()
k = 0
for (i in c(1:length(repl_pattern))) { # change after test-phase !!!
  tech_reps <- as.vector(unlist(repl_pattern[i]))
  sample_name <- names(repl_pattern)[i]
  for (j in 1:length(tech_reps)) {
    k = k + 1
    repl1_nr <- read.table(TIC_files[j])
    bad_color_pos <- tech_reps[j] %in% remove_pos[[1]]
    bad_color_neg <- tech_reps[j] %in% remove_neg[[1]]
    if (bad_color_neg & bad_color_pos) {
	plotcolor = '#F8766D'
    } else if (bad_color_pos) {
	plotcolor = "#ED8141"
    } else if (bad_color_neg) {
	plotcolor = "#BF80FF"
    } else {
	plotcolor = 'white'
    }
    tic_plot <- ggplot(repl1_nr, aes(retentionTime, TIC)) +
      geom_line(linewidth = 0.3) +
      geom_hline(yintercept = highest_tic_max, col = "grey", linetype = 2, linewidth = 0.3) +
      labs(x = 't (s)', y = 'TIC', title = paste0(tech_reps[j], "  ||  ", sample_name)) +
      theme(plot.background = element_rect(fill = plotcolor), axis.text = element_text(size = 4), axis.title = element_text(size = 4), plot.title = element_text(size = 6))
    tic_plot_list[[k]] <- tic_plot
  }

}
# create a layout matrix dependent on numer of replicates
layout <- matrix(1:(10 * nr_replicates), 10, nr_replicates, TRUE)

tic_plot_pdf <- marrangeGrob(grobs = tic_plot_list,
			     nrow = 10,
			     ncol = nr_replicates,
			     layout_matrix = layout,
			     top = quote(paste("TICs of run", run_name," \n colors: red = both modes misinjection, orange = pos mode misinjection, purple = neg mode misinjection \n ", g, "/", npages)))
ggsave(filename = paste0("./../../../", run_name, "_TICplots.pdf"), tic_plot_pdf, width = 21, height = 29.7, units = "cm")
