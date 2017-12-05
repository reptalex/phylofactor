#' Gets Groups considered at a given factor level
#' 
#' @export
#' @param PF phylofactor class object from \code{\link{PhyloFactor}}
#' @param factor factor level of interest, must be greater than 1. For factor=1, use \code{\link{getPhyloGroups}}. A list of paired groups considered at that level will be returned.
#' @return output similar to \code{\link{getPhyloGroups}}
pf.getPhyloGroups <- function(PF,factor){
  Grps <- phylofactor::getGroups(tree)
  PFgrps <- PF$groups[1:(factor-1)]
  tree <- PF$tree
  treeList <- list(tree)
  binList <- list(1:ape::Ntip(tree))
  
  for (grp in PFgrps){
    grp <- getLabelledGrp(tree=tree,Groups=grp)
    treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
    binList <- updateBinList(binList,grp)
    Grps <- getNewGroups(tree,treeList,binList)
  }
  return(Grps)
}