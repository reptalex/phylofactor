#'Input OTUs and outputs their taxonomic detail, optinally trimmed to their common taxonomic levels.
#'
#'@param otus Species IDs found in first column of Taxonomy
#'@param Taxonomy Taxonomy whose first column has OTUs/species IDs and whose second column has greengenes-compatible taxonomic strings
#'@param common.name Logical indicating whether or not to trim the output to the longest common prefix in the taxonomic strings of all otus

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
