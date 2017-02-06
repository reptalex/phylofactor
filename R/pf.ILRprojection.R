#' Projects Data onto basis defined by phylofactor object
#'@export
#' @param PF Phylofactor object. See \code{\link{PhyloFactor}}
#' @param Data Optional new data matrix, Must have rownames matching those in PF$Data
#' @param nfactors Number of factors used for projection. Default nfactors=2 for ordination-visualization

pf.ILRprojection <- function(PF,Data=NULL,nfactors=2){


  if (is.null(Data)){
    output <- lapply(PF$groups[1:nfactors],amalg.ILR,LogData=log(PF$Data)) %>% unlist %>% matrix(.,nrow=nfactors,byrow=T)
    colnames(output) <- colnames(PF$Data)
  }
  else {
    if(!all(rownames(Data) %in% rownames(PF$Data))){stop('Not all rownames in data are found in PF$tree')}
    Data <- Data[match(rownames(PF$Data),rownames(Data)),]
    output <- lapply(PF$groups[1:nfactors],amalg.ILR,LogData=log(Data)) %>% unlist %>% matrix(.,nrow=nfactors,byrow=T)
    colnames(output) <- colnames(Data)
    rownames(output) <- rownames(Data)
  }

  return(output)
}
