#' Produces cluster object with necessary libraries to perform phylofactorization
#' @export
#' @param ncores number of cores
#' @param ... optional input arguments for \code{\link{makeCluster}}

phyloFcluster <- function(ncores=2,...){
  #creates cluster with all the necessary packages and functions to use these functions in parallel
  cl <- parallel::makeCluster(ncores,...)
  parallel::clusterEvalQ(cl,{library(magrittr)
                             library(phylofactor)
                             library(biglm)})
  return(cl)
}
