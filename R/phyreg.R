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
  #internal function for phyloregPar, using bigglm for memory-efficiency
  #input list of Groups, will output Y,GLMs, and, if choice=='var', residual variance.
  Y <- lapply(X=Grps,FUN=amalgamate,Data=Data,method)
  GLMs <- lapply(X=Y,FUN = pglm,x=XX,frmla=frmla,...)
  Yhat <-  lapply(GLMs,predict,newdata=NULL)

  reg <- NULL
  if(choice=='var'){
    # We need to predict the data matrix & calculate the residual variance.
    reg$residualvar <- mapply(PredictAmalgam,Yhat,Grps,n,method,Pbasis=Pbasis,SIMPLIFY=F) %>%
                           sapply(.,FUN = residualVar,Data=Data)
  }

  reg$Y <- Y
  # reg$GLMs <- GLMs #the GLMs from lapply(X=Y,FUN = phyloreg,x=XX,frmla=frmla,...) take up a lot of memory and have been removed.


  reg$stats <- mapply(GLMs,FUN=getStats,Y) %>%
                  unlist %>%
                  matrix(.,,ncol=2,byrow=T) #contains Pvalues and F statistics
    rownames(reg$stats) <- names(Y)
    colnames(reg$stats) <- c('Pval','F')
  reg$Yhat <- Yhat

  gc()

  return(reg)
}
