#' Color bins on a plot of a phylogeny
#'
#' @export
#' @param X Either a phylofactor object or a list of bins. See \code{\link{bins}} or \code{\link{PhyloFactor}}
#' @param factor Integer. Bins will be the set of bins defined by factors up to the input \code{factor} in phylofactor object.
#' @param tree phylo-class tree. If input \code{X} is a list of bins, tree must be input, otherwise it will be extracted from phylofactor object.
#' @param n Size of total set to use for colorfcn - this enables multiple plots with differnent numbers of bins to have the same color order.
#' @param colorfcn Color function. Default is rainbow.
#' @param type type of phylogeny for plot.phylo. Default is unrooted.
#' @param show.tip.label Logical. Whether or not to show tip labels. Default is False.
#' @param legend.return Logical indicating whether or not to return legend labels and colors for making a legend.
#' @return plot.phylo where unique edges for each bin are different colors, with optional return of legend information
#' @examples
#' library(phylofactor)
#' data('FTmicrobiome')
#'
#' binPhyloPlot(FTmicrobiome$PF,factor=3)
binPhyloPlot <- function(X,factor,tree=NULL,n=NULL,colorfcn=rainbow,type='unrooted',show.tip.label=F,legend.return=F,...){

  if (class(X)=='phylofactor'){
    tree <- X$tree
    X <- bins(X$basis[,1:factor,drop=F])
  } else {
    if (is.null(tree)){stop('if input X is not a phylofactor object, must also input a phylo-class object in tree')}
  }

  otus <- sapply(X,FUN=function(ix,PF) tree$tip.label[ix],PF=PF)
  if (is.null(n)){n=length(X)}
  Cols <- colorfcn(n)
  EdgeCols <- rep('black',Nedge(tree))


  Edgs <- lapply(otus,FUN=function(x,tree) extractEdges(tree,x,type=3),tree=tree)

  for (nn in 1:length(X)){
    EdgeCols[Edgs[[nn]]] <- Cols[nn]
  }

  plot.phylo(tree,edge.color=EdgeCols,type=type,show.tip.label = show.tip.label,...)

  l <- list(X,sapply(as.list(1:length(X)),FUN=function(x) paste('bin',x)),Cols)
  names(l) <- c('bins','Legend','Colors')
  if (legend.return){
    return(l)
  }

}
