#' Returns the taxonomic IDs for all the OTUs in the atoms formed by splits in a Phylofactor object.
#' @param PF Phylofactor object
#' @param Taxonomy taxonomy list. First column must be species IDs corresponding to tip-labels of the tree, second column must be taxonomy strings (currently, this works for greengenes taxonomy strings - beware other taxonomy formats)
#' @param tree Phylogeny
#' @param clades Clades, or factors, of phylofactor object to be inlcuded. Currently, this must be a sequence from 1:n, i.e. cannot be c(2,3,4...,n) or have any missing clades. The reason for this is the atom function currently requires input of a sequential partition.
#' @param common.name Logical indicating whether or not to truncate the output to only one taxonomic assignment per atom, based on the finest, common name to all taxa in the atom.


phylofactor.TaxaIDs <- function(PF,Taxonomy,tree,clades=NULL,common.name=F){
  if (is.null(clades)){clades <- 1:length(PF$nodes)}
  #what are the atoms for a lower-dimensional factorization?
  atms <- atoms(PF$basis[,clades])
  #let's get the taxonomic IDs of all those atms
  atom.otus <- lapply(atms,FUN=function(x,tree){return(tree$tip.label[x])},tree)
  taxatoms <- lapply(atom.otus,FUN=OTUtoTaxa,Taxonomy=Taxonomy,common.name=common.name)
  return(taxatoms)
}
