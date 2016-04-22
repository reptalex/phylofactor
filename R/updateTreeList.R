#' internal Phylofactor function - splits the tree corresponding to the split atom
#'
#' @export
#' @param treeList list of trees
#' @param atomList list of atoms corresponding to treeList
#' @param grp two-element list corresponding to the splitting of an atom
#' @examples
#' library(ape)
#' set.seed(6)
#' tree <- rtree(10)
#' Groups <- getGroups(tree)
#' treeList <- list(tree)
#' atomList <- list(1:10)
#' factor <- 16
#' grp <- getLabelledGrp(factor,tree,Groups)
#' grp
#' lapply(grp,FUN=function(g,tree){tree$tip.label[g]},tree=tree)
#'
#' node=17
#'
#' treeList <- updateTreeList(treeList,atomList,grp,tree)
#' par(mfrow=c(1,3))
#' plot.phylo(tree,main='Original Tree',cex=2)
#' edgelabels('SPLIT',12,bg='red',cex=2)
#'
#' plot.phylo(treeList[[1]],main='First Atom Tree',cex=2)
#' plot.phylo(treeList[[2]],main='Second Atom Tree',cex=2)

updateTreeList <- function(treeList,atomList,grp,tree){
  ix <- whichAtomSplit(grp,atomList)
  tips=sum(grepl('tip',names(grp))) ## How many tips are there? They will be removed from AtomList and recalculated later.
  tr <- treeList[[ix]]  #this is the sub-tree being split
  leaves <- tr$tip.label

  if (all.equal(tr,tree)==FALSE){
    tmap <- match(tr$tip.label,tree$tip.label)
    grp <- lapply(grp,FUN = function(grp,tmap) match(grp,tmap),tmap=tmap)
  }

  if (tips==0){
    treeList[[ix]] <- drop.tip(tr,setdiff(leaves,leaves[grp[[1]]]))
    treeList[[length(treeList)+1]] <- drop.tip(tr,setdiff(leaves,leaves[grp[[2]]]))
  } else if (tips==1){
    cld <- setdiff(1:2,which(grepl('tip',names(grp))))
    #replace the split tree with the clade
    treeList[[ix]] <- drop.tip(tr,setdiff(leaves,leaves[grp[[cld]]]))
  } else {
    treeList[[ix]] <- NULL
  }
  return(treeList)
}
