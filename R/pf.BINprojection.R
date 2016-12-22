#' Projects phylofactored data onto bins defined by the factors 1:factor
#' 
#' @export
#' @param PF phylofactor object. See \link{PhyloFactor}
#' @param factor integer between 0 and \code{PF$nfactors}
#' @param taxonomy Optional taxonomy. If input, the list of OTUs in each bin will be the taxonomic names, not the OTU IDs.
#' @param common.name Logical. If input taxonomy, will indicate whether to trim taxonomic list to the longest common prefix of each bin.
#' @param uniques Logical. If input taxonomy, will indicate whether to trim the taxonomic list to the unique names.
#' @param prediction Logical. If true, will output the predictions for each bin instead of the observations.
#' @param rel.abund Logical. If true, will return the total relative abundances in the bin. Default, F, returns the geometric mean of the relative abundances of each bin
#' @return Returns list containing the compositional dataset formed by the bins and the list of OTUs in each bin.


pf.BINprojection <- function(PF,factor=PF$nfactors,taxonomy=NULL,common.name=F,uniques=T,prediction=F,rel.abund=F){
  
  Bins <- bins(PF$basis[,1:factor])
  if (prediction){
    if (rel.abund){
      binned_Data <- lapply(Bins,FUN=function(ix,Y) colSums(Y[ix,,drop=F]),Y=pf.predict(PF,factors=factor))
    } else {
      binned_Data <- lapply(Bins,FUN=function(ix,Y) compositions::geometricmeanCol(Y[ix,,drop=F]),Y=pf.predict(PF,factors=factor))
    }
  } else {
    if (rel.abund){
      binned_Data <- lapply(Bins,FUN=function(ix,Y) colSums(Y[ix,,drop=F]),Y=PF$Data)
    } else {
      binned_Data <- lapply(Bins,FUN=function(ix,Y) compositions::geometricmeanCol(Y[ix,,drop=F]),Y=PF$Data)
    }
  }
  output <- NULL
  output$Data <- matrix(unlist(binned_Data),nrow=factor+1,byrow=T)
  output$Data <- t(t(output$Data)/colSums(output$Data))
  output$otus <- lapply(Bins,FUN=function(ix,Data) rownames(Data[ix,,drop=F]),Data=PF$Data)
  names(output$otus) <- sapply(as.list(1:(factor+1)),FUN = function(x) paste('Bin',x))
  if (!is.null(taxonomy)){
    output$otus <- lapply(output$otus,FUN = function(otu,taxonomy,common.name,uniques) OTUtoTaxa(otu,taxonomy,common.name,uniques),taxonomy=taxonomy,common.name=common.name,uniques=uniques)
  }
  rownames(output$Data) <- names(output$otus)
  colnames(output$Data) <- colnames(PF$Data)
  return(output)
}