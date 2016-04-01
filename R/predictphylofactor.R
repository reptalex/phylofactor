predict.phylofactor <- function(Factor,factors=NULL,...){
  #outputs the predictions, using nfactors from PhyloFactor.
  #additional arguments '...' are for function predict()

  if (is.null(factors)){factors=1:length(Factor$nodes)}

    coefs <- NULL
    d=1
    for (nn in factors){
      gg <- Factor$glms[[nn]]
      pp <- predict(gg,...)
      coefs <- matrix(c(coefs,pp),ncol=d,byrow=F)
      d=d+1
    }

    Dat <- coefs %>% apply(MARGIN=1,ilrInv,V=Factor$basis[,factors])

  return(Dat)
}
