#' internal function for parallelized phylofactor
#' 
#' @param nset set of nodes to be evaluated by worker
#' @param tree_map mapping of nodes to \code{treeList}
#' @param treetips number of tips in each tree
#' @param treeList list of phylo class objects.
getAssignedGroups <- function(nset,tree_map,treetips,treeList){
  Grps <- vector(mode='list',length=0)
  ix = 0
  for (nn in nset){
    ############# getting grp - list of group & complement - for node in nset #####
    if (nn>tree_map[1]){
      whichTree <- max(which(tree_map<nn))+1
      if ((nn-tree_map[[whichTree-1]])==treetips[whichTree]+1){ 
        #This prevents us from drawing the root of a subtree, which has no meaningful ILR transform. 
        next 
      }
      grp[[1]] <- phangorn::Descendants(treeList[[whichTree]],node=(nn-tree_map[whichTree-1]))[[1]]
    } else {
      whichTree <- 1
      if (nn==(treetips[1]+1)){ 
        #This prevents us from drawing the root of a subtree, which has no meaningful ILR transform. 
        next 
      }
      grp[[1]] <- phangorn::Descendants(treeList[[1]],node=nn)[[1]]
    }
    grp[[2]] <- setdiff(1:treetips[whichTree],grp[[1]])
    ix = ix+1
    Grps[[ix]] <- lapply(grp,FUN=function(x,tree) tree$tip.label[x],tree=treeList[[whichTree]])
    #Converting numbered grps to species ensures we don't get mixed up.
  }
}