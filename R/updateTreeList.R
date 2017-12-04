#' internal Phylofactor function - splits the tree corresponding to the split bin
#'
#' @export
#' @param treeList list of trees
#' @param binList list of bins corresponding to treeList
#' @param grp two-element list corresponding to the splitting of an bin
#' @param tree phylo class object
#' @param skip.check Logical whether to skip checking if we grabbed the correct tree.
#' @examples
#' library(ape)
#' set.seed(6)
#' tree <- rtree(10)
#' Groups <- getGroups(tree)
#' treeList <- list(tree)
#' binList <- list(1:10)
#' factor <- 2
#' grp <- getLabelledGrp(factor,tree,Groups)
#' grp
#' lapply(grp,FUN=function(g,tree){tree$tip.label[g]},tree=tree)
#'
#' node=17
#'
#' treeList <- updateTreeList(treeList,binList,grp,tree)
#' par(mfrow=c(1,3))
#' plot.phylo(tree,main='Original Tree',cex=2)
#' edgelabels('SPLIT',factor,bg='red',cex=2)
#'
#' plot.phylo(treeList[[1]],main='First Bin Tree',cex=2)
#' plot.phylo(treeList[[2]],main='Second Bin Tree',cex=2)

updateTreeList <- function(treeList,binList,grp,tree,skip.check=F){
  ix <- phylofactor:::whichBinSplit(grp,binList)
  tips=sum(grepl('tip',names(grp))) ## How many tips are there? They will be removed from BinList and recalculated later.
  tr <- treeList[[ix]]  #this is the sub-tree being split
  leaves <- tr$tip.label

  if (skip.check){
    tmap <- match(tr$tip.label,tree$tip.label)
    grp <- lapply(grp,FUN = function(grp,tmap) match(grp,tmap),tmap=tmap)
  } else {
    if (all.equal(tr,tree)==FALSE){
    tmap <- match(tr$tip.label,tree$tip.label)
    grp <- lapply(grp,FUN = function(grp,tmap) match(grp,tmap),tmap=tmap)
    }
  }

  if (tips==0){
    treeList[[ix]] <- ape::drop.tip(tr,setdiff(leaves,leaves[grp[[1]]]))
    treeList[[length(treeList)+1]] <- ape::drop.tip(tr,setdiff(leaves,leaves[grp[[2]]]))
  } else if (tips==1){
    cld <- setdiff(1:2,which(grepl('tip',names(grp))))
    #replace the split tree with the clade
    treeList[[ix]] <- ape::drop.tip(tr,setdiff(leaves,leaves[grp[[cld]]]))
  } else {
    treeList[[ix]] <- NULL
  }
  return(treeList)
}
