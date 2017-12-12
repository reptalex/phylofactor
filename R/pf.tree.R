#' ggtree-based plotting of phylofactor bins or factors
#'
#' @export
#' @param pf phylofactor class object
#' @param method either "factors" or "bins", based on whether to use factors ( \code{factors} along with \code{groups}) or resultant bins (\code{pf$bins}) from phylofactorization.
#' @param factors vector of factor numbers to be plotted
#' @param groups vector of groups (1 or 2) to be plotted for each factor. Default is 1
#' @param colors colors for each \code{factors x groups} pair, or each member of \code{GroupList} 
#' @param GroupList optional list of groups for plotting
#' @param bg.color background color
#' @param bg.alpha background alpha
#' @param alphas vector, of length = \code{nrow(factor.map)} or \code{length(GroupList)} or, if default settings, of length \code{pf$nfactors}. Controls alphas for highlights of each clade. Warning: the current output legend gives you the colors layered onto each clade. Layering transparent colors on top of one-another will produce a different color not in the legend..
#' @param layout See \code{\link{ggtree}}. Default is "circular"
#' @param rootnode Logical. If true, will fill in the root node. Default controls root only with bg.color 
#' @param color.fcn Color function for plotting. Default is rainbow.
#' @param ... optional arguments for ggtree
#' @examples 
#' library(phylofactor)
#' data(FTmicrobiome)
#' pf <- FTmicrobiome$PF
#' gg <- pf.tree(pf,factors=1:3,layout='circular')
#' gg$ggplot
#' gg$legend
pf.tree <- function(pf,method='factors',factors=NULL,groups=NULL,colors=NULL,GroupList=NULL,bg.color=NA,bg.alpha=0.1,alphas=NULL,layout='circular',rootnode=FALSE,top.layer=F,top.alpha=0.1,color.fcn=viridis::viridis,...){
  if (!(is.null(GroupList) & is.null(factors))){
    if (method=='bins'){
      warning('input GroupList or factors will override method="bins"')
      method='factors'
    }
  }
  if (method=='factors'){
    if (is.null(GroupList)){
      if (is.null(factors)){
        if (is.null(groups)){
          factor.map=data.frame('factor'=1:pf$nfactors,'group'=rep(1,pf$nfactors))
        } else {
          factor.map=data.frame('factor'=1:pf$nfactors,'group'=groups)
        }
      } else {
        if (is.null(groups)){
          factor.map <- data.frame('factor'=factors,'group'=1)
        } else {
          if (length(groups) != length(factors)){
            stop('length of groups is not equal to length of factors')
          } else {
            factor.map <- data.frame('factor'=factors,'group'=groups)
          }
        }
      }
      m <- nrow(factor.map)
      if (is.null(colors)){
        cols <- color.fcn(m)
        factor.map$colors <- cols
      } else {
        factor.map$colors <- colors
      }
    } else {
      m <- length(GroupList)
      method='groups'
      if (is.null(colors)){
        cols <- color.fcn(m)
      } else {
        cols <- colors
      }
      factor.map <- data.frame('Groups'=1:m,'colors'=cols)
    }
  }
  
  if (method=='bins'){
    GroupList <- pf$bins
    method='groups'
  }
  if (is.null(alphas)){
    alphas <- rep(1,max((pf$nfactors+1),length(GroupList)))
  }
  
  n=ape::Ntip(pf$tree)
  
  if (is.null(names)){names <- sapply(1:m,FUN=function(s) paste('Clade',s))}
  nd <- numeric(m)
  for (i in 1:(m)){
    if (method=='factors'){
      grp <- pf$groups[[factor.map[i,1]]][[factor.map[i,2]]]
      if (length(grp)>1){
        nd[i] <- ggtree::MRCA(pf$tree,pf$tree$tip.label[grp])
      } else {
        nd[i] <- grp
      }
    } else {
      grp <- GroupList[[i]]
      if (length(grp)>1){
        nd[i] <- ggtree::MRCA(pf$tree,pf$tree$tip.label[grp])
      } else {
        nd[i] <- grp
      }
    }
  }
  
  
  ix <- order(nd)                ## this gives us a map of which node to which factor/bin
  nd <- sort(nd,decreasing = F)
  
  gg <- ggtree::ggtree(pf$tree,layout=layout,...)
  if (nd[1]==(n+1) | !is.na(bg.color)){
    if (is.na(bg.color)){
      bg.color=cols[1]
    }
    gg <- gg+ggtree::geom_hilight(n+1,fill=bg.color,alpha=bg.alpha)
  }
  
  i=0
  for (ndd in nd){
    i=i+1
    if (!ndd==(n+1) | rootnode){
      gg <- gg+ggtree::geom_hilight(ndd,fill=cols[ix[i]],alpha=alphas[i])
    }
  }
  
  gg <- gg+ggtree::geom_tree(...)
  
  if (top.layer){
    i=0
    for (ndd in nd){
      i=i+1
      if (!ndd==(n+1) | rootnode){
        gg <- gg+ggtree::geom_hilight(ndd,fill=cols[ix[i]],alpha=top.alpha)
      }
    }
  }

  return(list('legend'=factor.map,'ggplot'=gg))
}