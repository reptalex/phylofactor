#' Outputs the atoms made by the first n factors in PhyloFactor object
#'
#'  @param PF Phylofactor object. See \code{\link{PhylOFactor}}
#'  @param nfactors Number of factors, from 1:nfactors, whose resultant atoms are desired.
#'  @return a list of atoms
atoms.PhyloFactor <- function(PF,nfactors=1){
  return(atoms(PF$ilrs[,1:nfactors]))
}
