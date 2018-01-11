#' Converts input data frame for \link{gpf} to matrix for faster indexing
#' 
#' @export
#' @param DF data frame containing, at the minimum, columns \code{"Species"}, \code{"Sample"}, and the column with character string defined in \code{mat.data}. If \code{DF} does not contain \code{"Sample"}, a vector will be output.
#' @param mat.data Character string. Column of \code{DF} input to the matrix
#' @param empty.val Value for empty Species x Sample pairs. Default is 0.
#' @examples 
#' library(phylofactor)
#' m=10   #number of species
#' n=15   #number of samples
#' d=100  #number of data points
#' DF <- data.frame('Species'=sample(letters[1:m],d,replace=T),
#'                  'Sample'=sample(1:n,d,replace=T),
#'                  'Abundance'=rpois(d,2))
#' M <- phyloframe.to.matrix(DF,mat.data='Abundance',empty.val=NA)
#' 
#' DF <- data.frame('Species'=letters[1:m],
#'                  'N'=1:m) 
#' phyloframe.to.matrix(DF)
phyloframe.to.matrix <- function(DF,mat.data='N',empty.val=0){
  samples <- unique(DF$Sample)
  if (is.null(samples)){
    samples <- 1
    DF$Samples <- 1
  }
  species <- unique(DF$Species)
  if (length(samples)==1 & any(table(DF$Species)>1)){
    stop('input DF does not have Samples, yet some species are replicated. Cannot assign multiple values to a species without there being multiple Samples.')
  }
  Mat <- matrix(empty.val,nrow=length(species),ncol=length(samples))
  rownames(Mat) <- as.character(species)
  colnames(Mat) <- as.character(samples)
  
  Mat[cbind(as.character(DF$Species),as.character(DF$Sample))] <- DF[[mat.data]]
  return(Mat)
}