#' internal generalized additive model method for \code{\link{PhyloFactor}} method='gam'
#' @param y ILR coordinate, y
#' @param X data frame, X
#' @param PF.output logical, if true will return customized output for phylofactor object
#' @param frmla formula for regression
#' @param choice Either "var" or "F". Determines which statistic to choose from gam for objective function
#' @param ... additional arguments for \code{\link{gam}}
#' @export
GAM <- function(y,X,PF.output=FALSE,frmla,choice,...){
  dataset <- cbind('Data'=y,X)
  gg <- mgcv::gam(frmla,data=dataset,...)
  
  if (PF.output){
    return(gg)
    break
  } else {
    output <- NULL
    if (choice=='var'){
      output$objective <- getStats(gg)['ExplainedVar']  
    } else {
      output$objective <- getStats(gg)['F']
    }
    output$stopStatistics <- getStats(gg)['Pval']
    return(output)
  }
}