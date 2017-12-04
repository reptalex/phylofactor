#' converts data frame in appropriate phylofactor format to matrix.
#' @export
#' @param DF data frame or data table containing columns "Sample", "Species", and "N". "N" is the observed data for generalized phylofactorization.
#' @examples 
#' library(phylofactor)
#' DF <- data.frame('Sample'=rep(1:3,each=3),'Species'=rep(c('a','b','c'),times=3),'N'=1:9)
#' phylo.frame.to.matrix(DF)
phylo.frame.to.matrix <- function(DF){
  samples <- unique(DF$Sample)
  if (is.null(samples)){
    samples <- 1
    DF$Samples <- 1
  }
  species <- unique(DF$Species)
  Mat <- matrix(0,nrow=length(species),ncol=length(samples))
  rownames(Mat) <- as.character(species)
  colnames(Mat) <- as.character(samples)
  
  Mat[cbind(as.character(DF$Species),as.character(DF$Sample))] <- DF$N
  return(Mat)
}
