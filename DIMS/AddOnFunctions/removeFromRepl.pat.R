removeFromRepl.pat <- function(bad_samples, repl_pattern, nr_replicates) {

  tmp = repl_pattern

  removeFromGroup=NULL

  for (i in 1:length(tmp)){
    tmp2 = repl_pattern[[i]]

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
    if (!is.null(remove)) repl_pattern[[i]]=repl_pattern[[i]][-remove]
  }

  if (length(removeFromGroup)!=0) {
    repl_pattern=repl_pattern[-removeFromGroup]
  }

  return(list("pattern"=repl_pattern))
}
