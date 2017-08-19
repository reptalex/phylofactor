#' ggtree-based plotting of phylofactor bins or factors
#'
#' @export
#' @param pf phylofactor class object
#' @param method either "factors" or "bins", based on whether to use factors (and \code{factor.map} or resultant bins from phylofactorization)
#' @param factor.map data frame with two columns. The first column is the factor number and the second is which group (1 or 2) should be highlighted
#' @param Grps optional list of groups for plotting
#' @param bg.color background color
#' @param bg.alpha background alpha
#' @param alphas vector, of length = \code{nrow(factor.map)} or \code{length(Grps)} or, if default settings, of length \code{pf$nfactors}. Controls alphas for highlights of each clade
#' @param layout. See \code{\link{ggtree}}. Default is "circular"
#' @param rootnode Logical. If true, will fill in the root node. Default controls root only with bg.color 
#' @param color.fcn Color function for plotting. Default is rainbow.
#' @param ... optional arguments for ggtree
#' @examples 
#' library(phylofactor)
#' data(FTmicrobiome)
#' PF <- FTmicrobiome$PF
#' gg <- pf.tree(PF,factor.map=data.frame('factors'=1:3,'groups'=rep(1,3)),layout='circular')
#' gg$ggplot
#' ## legend: PF$groups[gg$map] <--> gg$colors
pf.tree <- function(pf,method='factors',factor.map=NULL,Grps=NULL,bg.color=NA,bg.alpha=0.1,alphas=NULL,layout='circular',rootnode=FALSE,color.fcn=rainbow,...){
  if (!(is.null(Grps) & is.null(factor.map))){
    if (method=='bins'){
      warning('input Grps or factor.map will override metohd="bins"')
      method='factors'
    }
  }
  if (method=='factors'){
    if (is.null(Grps)){
      if (is.null(factor.map)){
        factor.map=data.frame('factor'=1:pf$nfactors,'group'=rep(1,pf$nfactors))
      }
      m <- nrow(factor.map)
    } else {
      m <- length(Grps)
      method='groups'
    }
  }
  if (method=='bins'){
    Grps <- pf$bins
    method='groups'
  }
  if (is.null(alphas)){
    alphas <- rep(0.5,max((pf$nfactors+1),length(Grps)))
  }
  
  n=Ntip(pf$tree)
  
  
  nd <- numeric(m)
  for (i in 1:(m)){
    if (method=='factors'){
      nd[i] <- MRCA(pf$tree,pf$tree$tip.label[pf$groups[[factor.map[i,1]]][[factor.map[i,2]]]])
    } else {
      nd[i] <- MRCA(pf$tree,pf$tree$tip.label[Grps[[i]]])
    }
  }
  cols <- color.fcn(m)
  
  ix <- order(nd)                ## this gives us a map of which node to which factor/bin
  nd <- sort(nd,decreasing = F)
  
  gg <- ggtree::ggtree(pf$tree,layout=layout,...)
  if (nd[1]==(n+1) | !is.na(bg.color)){
    gg <- gg+ggtree::geom_hilight(n+1,fill=bg.color,alpha=bg.alpha)
  }
  
  i=0
  for (ndd in nd){
    i=i+1
    if (!ndd==(n+1) | rootnode){
      gg <- gg+ggtree::geom_hilight(ndd,fill=cols[ix[i]],alpha=alphas[i])
    }
  }
  
  return(list('map'=ix,'colors'=cols,'ggplot'=gg))
}