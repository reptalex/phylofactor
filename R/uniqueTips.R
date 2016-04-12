#' Find all unique tips in a set of clades
#'
#' @param tree ape class phylogeny
#' @param Clades clades of interest
#' @return a vector of all the tips that are unique to \code{Clades} in \code{tree}
#' @examples
#' set.seed(1)
#' tree <- rtree(10)
#' plot.phylo(tree)
#' nodelabels()
#' Clades <- c(16,13)
#'
#' tips <- uniqueTips(tree,Clades)
#' tiplabels('***',tips,cex=3)


uniqueTips <- function(tree,Clades){
  Clds <- length(Clades)
  tips <- NULL
  for (cc in 1:Clds){
    tips <- sort(unique(c(tips,Descendants(tree,Clades[cc],type='tips')[[1]])))
  }
  return(tips)
}
