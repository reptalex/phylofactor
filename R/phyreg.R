#' Internal function for phyloregPar
#' @param Grps list of groups. see \code{\link{getGroups}}
#' @param Data data matrix. Rows are parts and columns are samples
#' @param XX independent variable
#' @param frmla formula for glm
#' @param n total number of taxa in Data
#' @param choice objective function allowing parallelization of residual variance calculations
#' @param ... optional inputs for glm

phyreg <- function(Grps,Data=NULL,Y=NULL,XX,frmla,n,choice,...){
  #internal function for phyloregPar, using bigglm for memory-efficiency
  #input list of Groups, will output Y,GLMs, and, if choice=='var', residual variance.
  reg <- NULL
  #Pre-allocation
  ngrps <- length(Grps)
  reg$Y <- vector(mode='list',length=ngrps)
  GLMs <- numeric(ngrps)
  reg$Yhat <- vector(mode='list',length=ngrps)
  
  if (is.null(Y)){
    reg$Y <- Grps %>% lapply(FUN=amalg.ILR,Log.Data=log(Data))
  } else {
    reg$Y <- Y
  }
  
  GLMs <- lapply(reg$Y,FUN = pglm,xx=XX,frmla=frmla,...)
  
  dataset <- c(list(rep(0,dim(Data)[2])),as.list(XX))
  names(dataset) <- c('Data',names(XX))
  dataset <- model.frame(frmla,data = dataset)
  reg$Yhat <- lapply(GLMs,predict,newdata=dataset)
  
  reg$stats <- mapply(GLMs,FUN=getStats,reg$Y) %>%
    unlist %>%
    matrix(.,ncol=3,byrow=T) #contains Pvalues and F statistics
  rownames(reg$stats) <- names(reg$Y)
  colnames(reg$stats) <- c('Pval','F','ExplainedVar')
  any.sigs=T
  if (!all(unlist(lapply(GLMs,FUN=function(bb){return(bb$converged)})))){warning('some bigGLMs did not converge. Try changing input argument maxit')}

  return(reg)
}
