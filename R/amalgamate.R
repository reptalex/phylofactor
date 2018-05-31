#' Amalgamation contrast function
#' 
#' @param Grps Either a list of \code{getPhyloGroups} or a single \code{getPhyloGroups} element
#' @param TransformedData Transformed Data. For amalgamation-phylofactorization, set \code{transform.fcn=I} and \code{contrast.fcn=amalgamate}.
#' @param pseudo.count Numeric greater than 0. If aggregation of a group is zero, it will be replaced with \code{pseudo.count}
#' @export
#' @examples 
#' library(phylofactor)
#' set.seed(1)
#' m=50
#' n=20
#' tree <- rtree(m)
#' Grps <- getPhyloGroups(tree)
#' M <- matrix(rpois(m*n,1),nrow=m)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- paste('Sample',1:n)
#' 
#' amalgamate(Grps,M)
#' amalgamate(Grps[[1]],M)
amalgamate <- function(Grps,TransformedData,pseudo.count=0.65){
  amlg <- function(grp,M)  log(sapply(colSums(M[grp[[1]],,drop=F]),max,pseudo.count)/sapply(colSums(M[grp[[2]],,drop=F]),max,pseudo.count))
  if (class(Grps[[1]])=='list'){
     return(t(sapply(Grps,amlg,TransformedData)))
  } else {
     return(amlg(Grps,TransformedData))
  }
}
