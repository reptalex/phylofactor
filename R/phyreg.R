#' Internal function for phyloregPar
#' @param Grps list of groups. see \code{\link{getGroups}}
#' @param Data data matrix. Rows are parts and columns are samples
#' @param XX independent variable
#' @param frmla formula for glm
#' @param n total number of taxa in Data
#' @param choice objective function allowing parallelization of residual variance calculations
#' @param method method for amalgamating groups
#' @param Pbasis coming soon
#' @param ... optional inputs for glm

phyreg <- function(Grps,Data,XX,frmla,n,choice,method,Pbasis,...){
  #internal function for phyloregPar
  #input list of Groups, will output Y,GLMs, and, if choice=='var', residual variance.
  Y <- lapply(X=Grps,FUN=amalgamate,Data=Data,method)
  GLMs <- lapply(X=Y,FUN = phyloreg,x=XX,frmla=frmla,...)
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
