#' replaces zeros and log-transforms data
logData <- function(Data){
  if (any(Data<0)){
    stop('Some Data<0. Cannot perform ilr transform on negative values.')
  }
  Data[Data==0] <- 0.65
  return(log(Data))
}