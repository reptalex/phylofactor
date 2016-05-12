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
#' 
#' 

TipTest <- function(PF,factors=1:(PF$nfactors),lower.tail=FALSE){
  nTips <- dim(PF$Data)[1]   ## number of tips
  nEdges <- 2*nTips-3               ## total number of edges
  nBasalEdges <- nEdges-nTips                 ## number of basal edges
  nFactors <- length(factors)
  expectedTips <- round(nFactors*(nTips/nEdges))
  obs <- sum(sapply(atoms(PF$basis[,factors,drop=F]),FUN=function(x) length(x)==1))
  
  obs <- sum(PF$factors[factors,1]=='tip' | PF$factors[factors,2]=='tip')
  
  if (nFactors==1){
    
    p <- obs*(nBasalEdges/nEdges)+(1-obs)*(nTips/nEdges)
    if (lower.tail){
      p <- 1-p
    }
    
  } else {
    p <- phyper(q=obs,m=nTips,n=nBasalEdges,k=nFactors,lower.tail=lower.tail)
  }
  
  
  ### Report the data
  tst <- c('upper tail hypergeometric test','lower tail hypergeometric test')[lower.tail+1]
  tbl <- matrix(c(nEdges,nTips,nBasalEdges,length(factors),expectedTips,obs,p),nrow=1)
  colnames(tbl) <- c('Nedges','Ntips','Nbasaledgs','Nfactors','expectedTips','ObsTips','P value')
  rownames(tbl) <- ''
  
  output <- list(tbl)
  names(output) <- tst
  return(output)
}
