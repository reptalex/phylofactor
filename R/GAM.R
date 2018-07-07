#' internal generalized additive model method for \code{\link{PhyloFactor}} method='gam'
#' @param y ILR coordinate, y
#' @param X data frame, X
#' @param PF.output logical, if true will return customized output for phylofactor object
#' @param gamfrmla formula for regression
#' @param gamchoice Either "var" or "F". Determines which statistic to choose from gam for objective function
#' @param ... additional arguments for \code{\link{gam}}
#' @export
GAM <- function(y,X,PF.output=FALSE,gamfrmla,gamchoice,...){
  dataset <- cbind('Data'=y,X)
  args=list('data'=dataset,'formula'=gamfrmla,...)
  gg <- do.call(mgcv::gam,args)
  
  if (PF.output){
    return(gg)
    break
  } else {
    output <- NULL
    if (gamchoice=='var'){
      output$objective <- getStats(gg)['ExplainedVar']  
    } else {
      output$objective <- getStats(gg)['F']
    }
    output$stopStatistics <- getStats(gg)['Pval']
    return(output)
  }
}