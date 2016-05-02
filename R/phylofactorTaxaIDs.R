#' Returns the taxonomic IDs for all the OTUs in the atoms formed by splits in a Phylofactor object.
#' @export
#' @param PF Phylofactor object
#' @param Taxonomy taxonomy list. First column must be species IDs corresponding to tip-labels of the tree, second column must be taxonomy strings (currently, this works for greengenes taxonomy strings - beware other taxonomy formats)
#' @param tree Phylogeny
#' @param factors factors of phylofactor object to be inlcuded. Currently, this must be a sequence from 1:n, i.e. cannot be c(2,3,4...,n) or have any missing factors. The reason for this is the atom function currently requires input of a sequential partition.
#' @param common.name Logical indicating whether or not to truncate the output to only one taxonomic assignment per atom, based on the finest, common name to all taxa in the atom.


phylofactor.TaxaIDs <- function(PF,Taxonomy,tree,nfactors=NULL,common.name=F){
  if (is.null(nfactors)){factors <- 1:PF$nfactors}
  #what are the atoms for a lower-dimensional factorization?
  atms <- atoms(PF$basis[,factors])
  #let's get the taxonomic IDs of all those atms
  atom.otus <- lapply(atms,FUN=function(x,tree){return(tree$tip.label[x])},tree)
  taxatoms <- lapply(atom.otus,FUN=OTUtoTaxa,Taxonomy=Taxonomy,common.name=common.name)

  names(taxatoms) <- lapply(atms,FUN=function(x,tree,nms) nms[1+as.numeric(ape::is.monophyletic(tree,x))],tree=tree,nms=c('Paraphyletic','Monophyletic'))


  return(taxatoms)
}
