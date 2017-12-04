#' replace zeros in Data for compositional data.
replaceZeros <- function(Data){
  Data[Data==0] <- 0.65
  return(Data)
}