#' Returns the taxonomic IDs for all the OTUs in the bins formed by splits in a Phylofactor object.
#' @export
#' @param PF Phylofactor object
#' @param Taxonomy taxonomy list. First column must be species IDs corresponding to tip-labels of the tree, second column must be taxonomy strings (currently, this works for greengenes taxonomy strings - beware other taxonomy formats)
#' @param tree Phylogeny
#' @param factors factors of phylofactor object to be inlcuded. Currently, this must be a sequence from 1:n, i.e. cannot be c(2,3,4...,n) or have any missing factors. The reason for this is the bin function currently requires input of a sequential partition.
#' @param common.name Logical indicating whether or not to truncate the output to only one taxonomic assignment per bin, based on the finest, common name to all taxa in the bin.


phylofactor.TaxaIDs <- function(PF,Taxonomy,tree,nfactors=NULL,common.name=F,uniques=F){
  if (is.null(nfactors)){factors <- 1:PF$nfactors}
  #what are the bins for a lower-dimensional factorization?
  atms <- bins(PF$basis[,1:nfactors])
  #let's get the taxonomic IDs of all those atms
  bin.otus <- lapply(atms,FUN=function(x,tree){return(tree$tip.label[x])},tree)
  taxbins <- lapply(bin.otus,FUN=OTUtoTaxa,Taxonomy=Taxonomy,common.name=common.name,uniques=uniques)

  names(taxbins) <- lapply(atms,FUN=function(x,tree,nms) nms[1+as.numeric(ape::is.monophyletic(tree,x))],tree=tree,nms=c('Paraphyletic','Monophyletic'))


  return(taxbins)
}
