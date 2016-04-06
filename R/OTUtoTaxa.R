OTUtoTaxa <- function(otus,Taxonomy,common.name=T){

    if(typeof(otus)!='character'){otus <- as.character(otus)}
      notus <- length(otus)
      ind <- match(otus,Taxonomy[,1])
      taxa <- sort(unlist(lapply(as.list(Taxonomy[ind,2]),FUN=toString)))

  if(common.name==T){
    return(substr(taxa[1], start = 1, stop = lcprefix(taxa[1], taxa[length(taxa)])))
  } else {
    return(taxa)
  }
}
