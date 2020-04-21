#' ggtree-based plotting of phylofactor bins or factors
#'
#' @export
#' @param pf phylofactor class object
#' @param tree An alternative tree on which to plot the phylogenetic factors \code{pf$tree}.
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
#' @param top.layer Logical. If true, will superimpose hilights on top of tree. If false, tree will layer on top of hilights.
#' @param top.alpha Numeric. Alpha value for \code{top.layer}
#' @param color.fcn Color function for plotting. Default is viridis.
#' @param ... optional arguments for ggtree
#' @examples 
#' library(phylofactor)
#' data(FTmicrobiome)
#' pf <- FTmicrobiome$PF
#' big.tree <- FTmicrobiome$tree
#' tax <- FTmicrobiome$taxonomy
#' gg <- pf.tree(pf,factors=1:3,layout='circular')
#' gg$ggplot
#' gg$legend
#' 
#' ## We can also focus on a universal tree
#' # pf.tree(pf,big.tree,factors=1:3)
#' 
#' ## and a sub-tree
#' bacteroidetes <- tax$OTU_ID[grepl('p__Bacteroidetes',tax$taxonomy)] %>%
#'  intersect(pf$tree$tip.label)
#' bacteroidetes.tree <- drop.tip(big.tree,setdiff(big.tree$tip.label,bacteroidetes))
#' pf.tree(pf,bacteroidetes.tree,factors=setdiff(1:pf$nfactors,41))
#' ## factor 41 contains a large, paraphyletic group that encompases all of the 
#' ## bacteroidetes - this will color our entire tree purple.
pf.tree <- function(pf,tree=NULL,method='factors',factors=NULL,ignore.tips=TRUE,groups=NULL,colors=NULL,GroupList=NULL,bg.color=NA,bg.alpha=0.1,alphas=NULL,layout='circular',rootnode=FALSE,top.layer=F,top.alpha=0.1,color.fcn=viridis::viridis,...){
  
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
        cols <- colors
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
      factor.map <- data.frame('Groups'=1:m,'colors'=cols,stringsAsFactors = F)
    }
  }
  
  if (method=='bins'){
    GroupList <- pf$bins
    method='groups'
  }
  if (is.null(alphas)){
    alphas <- rep(1,max((pf$nfactors+1),length(GroupList)))
  }
  
  if (is.null(tree)){
    tree=pf$tree
  }
  
  n=ape::Ntip(tree)
  
  
  nd <- NULL
  for (i in 1:m){
    if (method=='factors'){
      grp <- pf$tree$tip.label[pf$groups[[factor.map[i,1]]][[factor.map[i,2]]]]
    } else {
      grp <- pf$tree$tip.label[GroupList[[i]]]
    }
    
    grp <- intersect(grp,tree$tip.label)
    if (length(grp)>0 & any(setdiff(tree$tip.label,grp) %in% pf$tree$tip.label)){
      if (length(grp)>1){
        nd <- c(nd,ggtree::MRCA(tree,grp))
      } else {
        nd <- c(nd,match(grp,tree$tip.label))
      }
    }
  }
  if (is.null(nd)){
    stop('None of the factors/groups could be mapped to the input tree')
  } 
  if (all(nd==(n+1))){
    stop('All factors/groups were mapped to the root of the input tree - no meaningful pf.tree to make')
  }
  
  ix <- order(nd,decreasing = F)                ## this gives us a map of which node to which factor/bin
  nd <- nd[ix]
  
  gg <- ggtree::ggtree(tree,layout=layout,...)
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
