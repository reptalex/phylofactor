#' Performs Phylogenetic Principal Components Analysis
#' 
#' @export
#' @param Data Positive-valued data matrix whose rownames are tip-labels in the input \code{tree}. 
#' @param tree phylo-class object whose tip-labels cover the rownames of Data.
#' @param ncores optional number of cores for built-in parallelization. Be cautious of memory - each worker is sent a copy of the dataset and tree, producing ncores+1 copies of the dataset & tree.
#' @param ncomponents integer indicating number of phylogenetic components to extract
#' @param output.edges Logical, whether or not to output the edges in the input \code{tree} corresponding to phylogenetic components
#' @param tol tolerance for compositional matrix. Rounding error in large datasets can lead to columns of compositional data not summing to 1. 
#' @param quiet Logical, whether or not to quiet warnings.
#' @param ... additional input arguments to \code{\link{PhyloFactor}}
#' @return PhyCA object containing Data, tree, basis, and edges of the phylogeny corresponding to each split.
#' @examples 
#' 
#' library(phylofactor)
#' data("FTmicrobiome")
#' 
#' Data <- FTmicrobiome$PF$Data
#' tree <- FTmicrobiome$PF$tree
#' X <- FTmicrobiome$X
#' taxonomy <- FTmicrobiome$taxonomy
#' clr <- function(A) apply(A,MARGIN=2,FUN=function(a) log(a)-mean(log(a)))
#' pf.heatmap(tree=tree,Data=clr(Data))
#' 
#' 
#' phca <- PhyCA(Data,tree,ncomponents = 2)
#' phcaPAR <- PhyCA(Data,tree,ncomponents=2,ncores=2)
#' 
#' pf.heatmap(tree=tree,Data=clr(Data))
PhyCA <- function(Data,tree,ncores=NULL,ncomponents=NULL,output.edges=T,tol=1e-5,quiet=T,...){
  
  
  output <- PhyloFactor(Data,tree,method='max.var',ncores=ncores,nfactors=ncomponents,tolerance = tol,quiet=quiet,...)
  
  if (output.edges){
    if (!is.null(ncores)){
      output$edges <- getFactoredEdgesPAR(PF=output,ncores=ncores)
      names(output$edges) <- sapply(as.list(1:length(output$groups)),FUN=function(a) paste('factor',a,sep=' '))
    } else {
      output$edges <- lapply(output$groups,FUN=function(grps,tree) getFactoredEdges(grp1=grps[[1]],grp2=grps[[2]],tree=tree),tree=output$tree)
    }
  }
  output$phylofactor.fcn <- 'PhyCA'
  class(output) <- 'phylofactor'
  return(output)
  
}
