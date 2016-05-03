


phylofactor.ATOMprojection <- function(PF,Data=NULL,nfactors=2,clr.transform=F){

  factors=1:nfactors
  atms <- atoms(-PF$basis[,1:nfactors,drop=F])

  if (is.null(Data)){
    output <- amalgamate(atms,PF$Data,method=PF$method,collapse=T)
  }
  else {
    if(!all(rownames(Data) %in% rownames(PF$Data))){stop('Not all rownames in data are found in Phylofactor tree')}
    dum <- PF$Data
    dum[match(rownames(Data),rownames(dum)),] <- Data
    output <- amalgamate(atms,PF$Data,method=PF$method,collapse=T)
    colnames(output) <- colnames(Data)
    rownames(output) <- rownames(Data)
  }
  if (clr.transform){
    output <- t(compositions::clr(t(output)))
  }
  return(output)
}
