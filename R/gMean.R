#' Computes geometric mean using logarithms
#'
#' @export
#' @param Y matrix or vector
#' @param MARGIN input for \code{\link{apply}}

gMean <- function(Y,MARGIN=2){
  g <- apply(Y,MARGIN=MARGIN,FUN=function(y) exp(mean(log(y))))
  return(g)
}
