phyloregPar <- function(Grps,Data,X,frmla,choice,method,Pbasis,cl,...){

n <- dim(Data)[1]
m <- length(Grps)
parG <- lapply(clusterSplit(cl,1:m),FUN <- function(ind,g){return(g[ind])},g=Grps)

reg <- parLapply(cl,X=parG,fun=phyreg,Data=Data,XX=X,frmla=frmla,n=n,choice=choice,method=method,Pbasis=Pbasis,...)

  output <- NULL
  for (pp in 1:length(parG)){
    output$GLMs <- c(output$GLMs,reg[[pp]]$GLMs)
    output$Y <- c(output$Y,reg[[pp]]$Y)
    output$residualvar <- c(output$residualvar,reg[[pp]]$residualvar)
  }
  
  return(output)
}