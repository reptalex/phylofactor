#' get bin number for species
#' 
#' @param B bin list. See \code{bins}
phylobin <- function(B){
  bn <- numeric(sum(sapply(B,length)))
  for (i in 1:length(B)){
    bn[B[[i]]] <- i
  }
  return(bn)
}