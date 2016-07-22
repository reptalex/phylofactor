#' Updates Bin list given a split group
#'
#' @export
#' @param binList List of bins
#' @param grp Group indicating the split of an bin
#' @examples
#' library(ape)
#' set.seed(6)
#' tree <- rtree(10)
#' Groups <- getGroups(tree)
#' treeList <- list(tree)
#' binList <- list(1:10)
#' factor <- 16
#' grp <- getLabelledGrp(factor,tree,Groups)
#' grp
#'
#' treeList <- updateTreeList(treeList,binList,grp,tree)
#' binList <- updateBinList(binList,grp)
#'
#'
#' Groups <- getNewGroups(tree,treeList,binList)
#'
#' factor=2
#' grp <- getLabelledGrp(factor,tree,Groups)
#' grp
#' unlist(grp)
#' tree$tip.label[unlist(grp)]
#'
#' treeList <- updateTreeList(treeList,binList,grp,tree)
#' binList <- updateBinList(binList,grp)
#'
#'
#'
#' Groups <- getNewGroups(tree,treeList,binList)
#' factor = 11
#' grp <- getLabelledGrp(factor,tree,Groups)
#' grp
#' unlist(grp)
#' tree$tip.label[unlist(grp)]
#'
#' treeList <- updateTreeList(treeList,binList,grp,tree)
#' binList <- updateBinList(binList,grp)
updateBinList <- function(binList,grp){
  ix <- whichBinSplit(grp,binList)
  tips=sum(grepl('tip',names(grp))) ## How many tips are there? They will be removed from BinList and recalculated later.

  if (tips==0){
    binList[[ix]] <- grp[[1]]
    binList[[length(binList)+1]] <- grp[[2]]
  } else if (tips==1){
    cld <- setdiff(1:2,which(grepl('tip',names(grp))))
    #replace the split bin with the clade
    binList[[ix]] <- grp[[cld]]
  } else {
    binList[[ix]] <- NULL
  }
  return(binList)
}
