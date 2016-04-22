#' Updates Atom list given a split group
#'
#' @export
#' @param atomList List of atoms
#' @param grp Group indicating the split of an atom
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
#'
#' treeList <- updateTreeList(treeList,atomList,grp,tree)
#' atomList <- updateAtomList(atomList,grp)
#'
#'
#' Groups <- getNewGroups(tree,treeList,atomList)
#'
#' factor=2
#' grp <- getLabelledGrp(factor,tree,Groups)
#' grp
#' unlist(grp)
#' tree$tip.label[unlist(grp)]
#'
#' treeList <- updateTreeList(treeList,atomList,grp,tree)
#' atomList <- updateAtomList(atomList,grp)
#'
#'
#'
#' Groups <- getNewGroups(tree,treeList,atomList)
#' factor = 11
#' grp <- getLabelledGrp(factor,tree,Groups)
#' grp
#' unlist(grp)
#' tree$tip.label[unlist(grp)]
#'
#' treeList <- updateTreeList(treeList,atomList,grp,tree)
#' atomList <- updateAtomList(atomList,grp)
updateAtomList <- function(atomList,grp){
  ix <- whichAtomSplit(grp,atomList)
  tips=sum(grepl('tip',names(grp))) ## How many tips are there? They will be removed from AtomList and recalculated later.

  if (tips==0){
    atomList[[ix]] <- grp[[1]]
    atomList[[length(atomList)+1]] <- grp[[2]]
  } else if (tips==1){
    cld <- setdiff(1:2,which(grepl('tip',names(grp))))
    #replace the split atom with the clade
    atomList[[ix]] <- grp[[cld]]
  } else {
    atomList[[ix]] <- NULL
  }
  return(atomList)
}
