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
  nn <- dim(PF$Data)[1]   ## number of tips
  N <- 2*nn-3               ## total number of edges
  mm <- N-nn                 ## number of basal edges
  nf <- length(factors)
  obs <- sum(PF$factors[factors,1]=='tip')

  if (nf==1){

      p <- obs*(mm/N)+(1-obs)*(nn/N)
      if (lower.tail){
        p <- 1-p
      }

  } else {
    p <- phyper(q=obs,m=nn,n=mm,k=nf,lower.tail=lower.tail)
  }

  tst <- c('upper tail hypergeometric test','lower tail hypergeometric test')[lower.tail+1]
  tbl <- matrix(c(N,nn,mm,length(factors),obs,p),nrow=1)
  colnames(tbl) <- c('Nedges','Ntips','Nbasaledgs','Nfactors','ObsTips','P value')
  rownames(tbl) <- ''

  output <- list(tbl)
  names(output) <- tst
  return(output)
}
