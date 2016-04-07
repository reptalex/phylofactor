phylofactor.TaxaIDs <- function(PF,Taxonomy,tree,clades=NULL,common.name=F){
  if (is.null(clades)){clades <- 1:length(PF$nodes)}
  #what are the atoms for a lower-dimensional factorization?
  atms <- atoms(PF$basis[,clades])
  #let's get the taxonomic IDs of all those atms
  atom.otus <- lapply(atms,FUN=function(x,tree){return(tree$tip.label[x])},tree)
  taxatoms <- lapply(atom.otus,FUN=OTUtoTaxa,Taxonomy=Taxonomy,common.name=common.name)
  return(taxatoms)
}
