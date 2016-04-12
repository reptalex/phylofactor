#' Projects Data onto basis defined by phylofactor object
#'
#' @param PF Phylofactor object. See \code{\link{PhyloFactor}}
#' @param Data Optional new data matrix, Must have rownames matching those in PF$Data
#' @param nfactors Number of factors used for projection. Default nfactors=2 for ordination-visualization

pf.projection <- function(PF,Data=NULL,nfactors=2){
  if (is.null(Data)){
    output <- t(ilr(t(PF$Data),V=PF$basis[,1:nfactors]))
    colnames(output) <- colnames(PF$Data)
  }
  else {
    if(!all(rownames(Data) %in% rownames(PF$Data))){stop('Not all rownames in data are found in Phylofactor tree')}
    dum <- PF$Data
    dum[match(rownames(Data),rownames(dum)),] <- Data
    output <- t(ilr(t(Data),V=PF$basis[1:nfactors]))
    colnames(output) <- colnames(Data)
    rownames(output) <- rownames(Data)
  }

  return(output)
}
