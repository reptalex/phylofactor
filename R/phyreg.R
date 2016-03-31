phyreg <- function(Grps,Data,XX,frmla,n,choice,method,Pbasis,...){
  #internal function for phyloregPar
  #input list of Groups, will output Y,GLMs, and, if choice=='var', residual variance.
  Y <- lapply(X=Grps,FUN=amalgamate,Data=Data,method)
  GLMs <- lapply(X=Y,FUN = phyloreg,x=XX,frmla=frmla,choice,...)
  Yhat <- lapply(GLMs,predict)
  
  reg <- NULL
  if(choice=='var'){
    predictions <- mapply(PredictAmalgam,Yhat,Grps,n,method,Pbasis=Pbasis,SIMPLIFY=F)
    residualvar <- sapply(predictions,FUN = residualVar,Data=Data)
    reg$residualvar <- residualvar
  }
  
  reg$Y <- Y
  reg$GLMs <- GLMs
  return(reg)
}
