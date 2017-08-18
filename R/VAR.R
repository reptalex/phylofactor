#' Computes variance for \code{\link{PhyloFactor}} method='max.var'
#' @param y ILR coordinate, y
#' @param X data frame, X
#' @param PF.output logical, if true will return customized output for phylofactor object
#' @export
VAR <- function(y,X,PF.output=FALSE){
  v <- var(y)
  if (PF.output){
    return(v)
    break
  } else {
    output <- NULL
    output$objective <- v  
    output$stopStatistics <- 0
    return(output)
  }
}