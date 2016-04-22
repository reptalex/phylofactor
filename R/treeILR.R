#' Computes ilr basis corresponding to a phylogeny.
#' 
#' @param phylo Ape class phylo object. Must be ultrametric - if not, we recommend trying phylo <- chronopl(phylo,lambda=1) before inputting.
#' @examples 
#' library(ape)
#' library(compositions)
#' tree <- rtree(10)
#' V <- tree.ilr
#' 
#' Data <- matrix(rpois(100,lambda=20),nrow=5)
#' Data[Data==0] <- 0.65
#' 
#' # Get the ILR coordinates for a dataset
#' X <- t(ilr(t(Data),V))
#' 
tree.ilr <- function(phylo){
  #input phylogeny. Output will be a list containing (1) V basis, (2) list of Groups, (3) mapping from groups to ilr
  V = gsi.buildilrBase(gsi.merge2signary(as.hclust(phylo)$merge))
  
  Groups <- clade.members.list(phylo)
  
  GroupVector <- function(node,tree){
    #create vector in which the two descendant clades of a given node are +/-1, and all others are zero.
    #furthermore, all members of the same sub-clade will have the same sign
    node=as.numeric(node)
    ntips <- length(tree$tip.label)
    GrpVec <- rep(0,ntips)
    desc.nodes <- Descendants(tree,node,type='children')
    for (i in 1:2){
      dum <- Descendants(tree,desc.nodes[i],type='tips')[[1]]
      GrpVec[dum]=sign(i-1.1)
    }
    return(GrpVec)
  }
  
  GrpVec <- lapply(X=as.list(names(Groups)),FUN = GroupVector,tree=phylo)
  
  VecMatch <- function(Vec,ILR){
    if (is.list(Vec)){
      Vec=Vec[[1]]
      warning('Coercing Vec out of list using [[1]]')
    }
    isalleq<-function(X,y){isTRUE(all.equal(X,y))}
    ix=which(apply(X=sign(ILR),MARGIN = 2,FUN = isalleq,y=Vec))
    if (length(ix)==0){
      Vec=-Vec
      ix=which(apply(X=sign(ILR),MARGIN = 2,FUN = isalleq,y=Vec))
    }
    return(ix)
  }
  
  Grp2ilr <- lapply(X=GrpVec,FUN=VecMatch,ILR=V)
  
  V <- V[,Grp2ilr]
  colnames(v) <- names(Groups)
  return(V)
}