#' Internal phylofactor function - get labelled Group from phylofactor
#' @export
#' @param factor integer indicating the element of \code{Groups} to be labelled
#' @param tree phylo object - global phylogeny used to indicate whether the elements of \code{Groups[factor]} are monophyletic
#' @param Groups list of group/complement pairs corresponding to atom splits
#' @examples
#' library(ape)
#' set.seed(6)
#' tree <- rtree(10)
#' Groups <- getGroups(tree)
#' treeList <- list(tree)
#' atomList <- list(1:10)
#' factor <- 16
#' grp <- getLabelledGrp(factor,tree,Groups)
#' grp

getLabelledGrp <- function(factor,tree,Groups){
  grp <- Groups[[factor]]

  # is it monophyletic?
  MonoPhy <- lapply(grp,FUN=function(tips,tree) as.integer(ape::is.monophyletic(tree,tips))+1,tree=tree)
  nms <- lapply(Type,FUN=function(a,nm) nm[a],nm=c('Paraphyletic','Monophyletic'))

  # is it a tip or a clade?
  Type <- lapply(grp,FUN=function(g) as.integer(length(g)==1)+1)
  typ <-  lapply(Type,FUN=function(a,nm) nm[a],nm=c('clade','tip'))

  names(grp) <- mapply(FUN=function(a,b) paste(a,b,sep=' '),a=nms,b=typ)
  if (any(names(grp)=='Monophyletic tip')){
    ix <- which(names(grp)=='Monophyletic tip')
    names(grp)[ix] <- 'tip'
  }
  return(grp)
}
