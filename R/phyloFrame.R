#' Tacks "phylo" variable onto data frame for \code{\link{gpf}}
#' 
#' @export
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
#' 
#' grp <- list(c(1:3),c(5:6))
#' phyloFrame(Data,grp,tree)
phyloFrame <- function(Data,grp,tree){
  factorFrame <- data.table('Species'=tree$tip.label,'phylo'='R')
  factorFrame[grp[[2]],'phylo'] <- 'S'
  ix <- setdiff(1:ape::Ntip(tree),unlist(grp))
  factorFrame[ix,'phylo'] <- NA
  data.table::setkey(factorFrame,Species)
  return(Data[factorFrame])
}
