#' internal version of getGroups to allow memory-efficient parallelization
#' 
#' @param tree phylo object
#' @param n node in \code{tree} whose descendants are used to get group/complement list.
#' @examples
#' 
#' library(ape)
#' set.seed(1)
#' tree <- rtree(15)
#' plot.phylo(tree,show.tip.label=F,use.edge.length=F)
#' tiplabels()
#' nodelabels()

getGroupsPar <- function(tree,n){
  output <- vector(mode='list',length=2)
  output[[1]] <- phangorn::Descendants(tree,n)[[1]]
  output[[2]] <- setdiff(1:(ape::Ntip(tree)),output[[1]])
  return(output)
}