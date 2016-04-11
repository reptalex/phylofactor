pf.projection <- function(PF,nfactors=2){

  output <- t(ilr(t(PF$Data),V=PF$basis[,1:nfactors]))
  colnames(output) <- colnames(PF$Data)
}
