#' replace zeros in Data for compositional data.
#' @param Data matrix or vector
replaceZeros <- function(Data){
  Data[Data==0] <- 0.65
  return(Data)
}