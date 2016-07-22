#' Color bins on a plot of a phylogeny
#'
#' @export
#' @param atms Either a list of bins or a phylofactor object. See \code{\link{bins}} or \code{\link{PhyloFactor}}
#' @param factor Integer. Bins will be the set of bins defined by factors up to the input \code{factor} in phylofactor object.
#' @param tree phylo-class tree. If input \code{atms} is a list of bins, tree must be input, otherwise it will be extracted from phylofactor object.
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
binPhyloPlot <- function(atms,factor,tree=NULL,n=NULL,colorfcn=rainbow,type='unrooted',show.tip.label=F,legend.return=F,...){

  if (class(atms)=='phylofactor'){
    tree <- atms$tree
    atms <- bins(atms$basis[,1:factor,drop=F])
  } else {
    if (is.null(tree)){stop('if input atms is not a phylofactor object, must also input a phylo-class object in tree')}
  }

  otus <- sapply(atms,FUN=function(ix,PF) tree$tip.label[ix],PF=PF)
  if (is.null(n)){n=length(atms)}
  Cols <- colorfcn(n)
  EdgeCols <- rep('black',Nedge(tree))


  Edgs <- lapply(otus,FUN=function(x,tree) extractEdges(tree,x,type=3),tree=tree)

  for (nn in 1:length(atms)){
    EdgeCols[Edgs[[nn]]] <- Cols[nn]
  }

  plot.phylo(tree,edge.color=EdgeCols,type=type,show.tip.label = show.tip.label,...)

  l <- list(atms,sapply(as.list(1:length(atms)),FUN=function(x) paste('bin',x)),Cols)
  names(l) <- c('bins','Legend','Colors')
  if (legend.return){
    return(l)
  }

}
