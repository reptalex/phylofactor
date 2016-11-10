#' returns the taxonomic identity of taxa in bins
#' @param PF Phylofactor object
#' @param nfactors number of factors used (will look at bins created by factors \code{1:nfactors})
#' @param taxonomy Taxonomy: first column is OTU IDs and second column is taxonomy.
#' @param common.name whether to simplify bin-taxonomy to its longest-similar-prefix
#' @param uniques - whether to include all unique taxa

binTaxa <- function(PF,nfactors,taxonomy,common.name=F,uniques=T){
  atms <- bins(PF$basis[,1:nfactors,drop=F])
  otus <- lapply(atms,FUN = function(x,names) names[x],names=PF$tree$tip.label)
  taxa <- lapply(otus,FUN=OTUtoTaxa,Taxonomy=taxonomy,common.name=common.name,uniques=uniques)
  taxa2 <- taxa
  
  if (uniques){
    for (nn in 1:(nfactors+1)){
      taxa2[[nn]] <- unique(uniqueTaxa(taxa[[nn]],unlist(taxa[setdiff(1:length(taxa),nn)])))
    }
    taxa=taxa2
  }
  
  return(taxa)
}