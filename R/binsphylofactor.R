#' Outputs the bins made by the first n factors in PhyloFactor object
#'  @param PF Phylofactor object. See \code{\link{PhylOFactor}}
#'  @param nfactors Number of factors, from 1:nfactors, whose resultant bins are desired.
#'  @return a list of bins
bins.PhyloFactor <- function(PF,nfactors=1){
  return(bins(PF$ilrs[,1:nfactors,drop=F]))
}
