#' Parallelized version of getFactoredEdges
#' 
#' @export
#' @param cl optional input \code{phyloFcluster}
#' @param ncores optional integer input ncores- must input wither cl or ncores.
#' @param tree phylo object
#' @param V partition matrix whose columns contain +/- 1 indicating group membership
#' @param PF phylofactor class object.
#' @examples 
#' data("FTmicrobiome")
#' V <- FTmicrobiome$PF$basis
#' tree <- FTmicrobiome$PF$tree
#' 
#' cl <- phyloFcluster(2)
#' 
#' getFactoredEdgesPAR(cl,tree=tree,V=V[,1:10])
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
  M <- ape::mrca(tree)
  
  splitFCN <- function(v,tree,M){
    apply(v,MARGIN=2,as.list)  %>% 
      lapply(.,FUN=function(v,tree,M) getFactoredEdges(v=v,tree=tree,M=M),tree=tree,M=M)
  }
  
  edgsPAR <- parallel::parLapply(cl,v,function(v,tree,M) splitFCN(v,tree,M),tree=tree,M=M)
  
  edgs <- vector(mode='list',length=ncol(V))
  for (pp in 1:length(ix)){
    edgs[ix[[pp]]] <- edgsPAR[[pp]]
  }
  
  if (!is.null(cl)){
    parallel::stopCluster(cl)
    rm('cl')
  }
  
  return(edgs)
}