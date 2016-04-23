#' Internal phylofactor function - get new grouplist corresponding to atomlist & treelist
#' @export
#' @param tree phylo object - global tree
#' @param treeList list of subtrees for each atom in \code{atomList}
#' @param atomList list of atoms corresponding to current factorization
#' @return List of groups for amalgamation. For example, see \code{updateAtomList}
getNewGroups <- function(tree,treeList,atomList){

  g = lapply(treeList,FUN=getGroups)
  tmap <- lapply(treeList,function(a,b) match(a$tip.label,b$tip.label),b=tree)
  #Now we need to write g - which has indexes corresponding to its treeList tip-labels - to have indexes corresponding to the global tree tiplabels
  g <- mapply(FUN=function(g,tmap) lapply(g,function(g,tmap) lapply(g,function(g,tmap) tmap[g],tmap=tmap),tmap=tmap),tmap=tmap,g=g)
  names(g) <- NULL
  if (is.null(dim(g))){
    G <- unlist(g,recursive=F)
  } else {
    m <- dim(g)[2]
    G <- NULL
    for (mm in 1:m){
      G <- c(G,g[,mm])
    }
  }
  ls <- sapply(G,length)
  G <- G[ls>0]
  return(G)
}
