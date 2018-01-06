#' get bin number for species
#' 
#' @param B bin list. See \code{bins}
#' @param pf phylofactor object - if input, will override input \code{B}
#' @param nfactors number of factors to be used - will result in \code{nfactors+1} bins
phylobin <- function(B,pf=NULL,nfactors=NULL){
  
  if (!is.null(pf)){
    if (is.null(nfactors)){
      nfactors=pf$nfactors
    } else {
      if (nfactors > pf$nfactors){
        stop('input nfactors is greater than pf$nfactors')
      }
    }
    B <- bins(pf$basis[,1:nfactors])
  }
  
  bn <- numeric(sum(sapply(B,length)))
  for (i in 1:length(B)){
    bn[B[[i]]] <- i
  }
  return(bn)
}