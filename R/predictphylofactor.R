predict.phylofactor <- function(Factor,factors=NULL,method='ILR',...){
  #outputs the predictions, using nfactors from PhyloFactor.
  #additional arguments '...' are for function predict()

  if (is.null(factors)){factors=1:length(Factor$nodes)}

  if (method=='ILR'){
    coefs <- NULL
    d=1
    for (nn in factors){
      gg <- Factor$glms[[nn]]
      coefs <- matrix(c(coefs,predict(gg,...)),ncol=d,byrow=F)
      d=d+1
    }

    Dat <- coefs %>% apply(MARGIN=1,ilrInv,V=Factor$basis[,factors])

  } else { ### Coming Soon - predict for log-ratio, with additive amalgamation

  }

  return(Dat)
}
