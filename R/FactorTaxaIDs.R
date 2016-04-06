FactorTaxaIDs <- function(Factor,clades,Taxonomy,tree,common.name=F){
  #what are the atoms for a lower-dimensional factorization?
  atms <- atoms(Factor$basis[,clades])
  #let's get the taxonomic IDs of all those atms
  atom.otus <- lapply(atms,FUN=function(x,tree){return(tree$tip.label[x])},tree)
  taxatoms <- lapply(atom.otus,FUN=OTUtoTaxa,Taxonomy=Taxonomy,common.name=common.name)
  return(taxatoms)
}
