atoms.PhyloFactor <- function(PF,nfactors=1){
  return(atoms(PF$ilrs[,1:nfactors]))
}
