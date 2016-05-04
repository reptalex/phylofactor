#' Test of whether tips are over/under-represented in phylofactorization
#' @export
#' @param PF phylofactor object. See \code{\link{PhyloFactor}}
#' @param factors Set of factors to test. Default is all the factors in PF input.
#' @param two.tailed Logical indicating if test is to be two-tailed.
#' @description Performs hypergeometric test on the number of tips in a phylofactorization
#' @examples
#' data('FTmicrobiome')
#'
#' TipTest(Ftmicrobiome$PF)

TipTest <- function(PF,factors=1:(PF$nfactors),lower.tail=FALSE){
  n <- dim(PF.F$Data)[1]   ## number of tips
  N <- 2*n-3               ## total number of edges
  m <- N-n                 ## number of basal edges
  nf <- length(factors)
  obs <- sum(PF.F$factors[factors,1]=='tip')

  p <- phyper(obs,n,m,length(factors),lower.tail=lower.tail)

  tst <- c('upper tail hypergeometric test','lower tail hypergeometric test')[lower.tail+1]
  tbl <- matrix(c(N,n,m,obs,length(factors),p),nrow=1)
  colnames(tbl) <- c('Nedges','Ntips','Nbasaledgs','Nfactors','ObsTips','P value')
  rownames(tbl) <- ''

  output <- list(tbl)
  names(output) <- tst
  return(output)
}
