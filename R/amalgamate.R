#' Amalgamation contrast function
#' 
#' @param Grps Either a list of \code{getPhyloGroups} or a single \code{getPhyloGroups} element
#' @param TransformedData Transformed Data. For amalgamation-phylofactorization, set \code{transform.fcn=I} and \code{contrast.fcn=amalgamate}.
#' @export
#' @examples 
#' library(phylofactor)
#' set.seed(1)
#' m=10
#' n=5
#' tree <- rtree(m)
#' Grps <- getPhyloGroups(tree)
#' M <- matrix(rlnorm(m*n),nrow=m)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- paste('Sample',1:n)
#' 
#' amalgamate(Grps,M)
#' amalgamate(Grps[[1]],M)
amalgamate <- function(Grps,TransformedData){
  amlg <- function(grp,M)  log(colSums(M[grp[[1]],,drop=F])/colSums(M[grp[[2]],,drop=F]))
  if (class(Grps[[1]])=='list'){
     return(t(sapply(Grps,amlg,TransformedData)))
  } else {
     return(amlg(Grps,TransformedData))
  }
}
