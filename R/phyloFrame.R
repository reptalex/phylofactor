#' Tacks "phylo" variable onto data frame for \code{\link{gpf}}
#' 
#' @param Data data table input to \code{\link{gpf}}. Must have column \code{"Species"}
#' @param grp Group list element from \code{\link{getPhyloGroups}}
#' @param tree phylo class object
#' @examples
#' library(phylofactor)
#' set.seed(1)
#' tree <- rtree(10)
#' Data <- data.table('Species'=tree$tip.label,'x'=rnorm(10),key='Species')
#' Grps <- getPhyloGroups(tree)
#' phylofactor:::phyloFrame(Data,Grps[[1]],tree)
phyloFrame <- function(Data,grp,tree){
  factorFrame <- data.table('Species'=tree$tip.label,'phylo'='R')
  factorFrame[grp[[2]],'phylo'] <- 'S'
  return(Data[factorFrame])
}
