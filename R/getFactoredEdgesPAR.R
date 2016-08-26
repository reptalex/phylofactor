#' Parallelized version of getFactoredEdges
#' 
#' @export
#' @param cl optional input \code{phyloFcluster}
#' @param ncores optional integer input ncores- must input wither cl or ncores.
#' @param tree phylo object
#' @param V partition matrix whose columns contain +/- 1 indicating group membership
#' @examples 
#' data("FTmicrobiome")
#' V <- FTmicrobiome$PF$basis
#' tree <- FTmicrobiome$PF$tree
#' 
#' cl <- phyloFcluster(2)
#' 
#' getFactoredEdges.par(cl,tree=tree,V=V[,1:10])
getFactoredEdgesPAR <- function(cl=NULL,ncores=NULL,tree=NULL,V=NULL,PF=NULL){
  if (is.null(cl)){
    if (is.null(ncores)){
      stop('must input either phyloFcluster cl or ncores to getFactoredEdgesPAR')
    } else {
      cl <- phylofactor::phyloFcluster(ncores)
    }
  }
  if (is.null(V)){
    if (is.null(PF)){
      stop('must input either partition matrix, V, or phylofactor object, PF')
    } else {
      if (class(PF)!='phylofactor'){
        stop('input PF must be of class phylofactor')
      }
      V <- PF$basis
      tree <- PF$tree
    }
  }
  
  if (is.null(tree)){
    stop('if inputting V, must input tree')
  }
  
  ind <- 1:ncol(V)
  ix <- parallel::clusterSplit(cl,ind)
  v <- lapply(ix,function(j,v) v[,j,drop=F],v=V)
  
  edgsPAR <- parallel::parLapply(cl,v,function(v,tree) lapply(apply(v,MARGIN=2,as.list),getFactoredEdges,tree=tree),tree=tree)
  
  edgs <- vector(mode='list',length=ncol(V))
  for (pp in 1:length(ix)){
    edgs[ix[[pp]]] <- edgsPAR[[pp]]
  }
  return(edgs)
}