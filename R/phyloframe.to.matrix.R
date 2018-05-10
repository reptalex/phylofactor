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
#' DF <- data.frame('Species'=sample(letters[1:m],d,replace=TRUE),
#'                  'Sample'=sample(1:n,d,replace=TRUE),
#'                  'Abundance'=rpois(d,2))
#' M <- phyloframe.to.matrix(DF,mat.data='Abundance',empty.val=0)
#' 
phyloframe.to.matrix <- function(DF,mat.data='N',empty.val=0){
  if (is.null(DF$Sample)){
    stop('must have element "Sample" for input DF')
  }
  if (is.null(DF$Species)){
    stop('must have element "Species" for input DF')
  }
  if (!all(mat.data %in% names(DF))){
    stop('not all mat.data are in names of DF')
  }
  
  samples <- unique(DF$Sample)
  species <- unique(DF$Species)
  Mat <- matrix(empty.val,nrow=length(species),ncol=length(samples))
  rownames(Mat) <- as.character(species)
  colnames(Mat) <- as.character(samples)
  
  if (length(mat.data)==1){
    Mat[cbind(as.character(DF$Species),as.character(DF$Sample))] <- DF[[mat.data]]
  } else {
    output <- vector(mode='list',length=length(mat.data))
    i=0
    for (mm in mat.data){
      i=i+1
      output[[i]] <- Mat
      output[[i]][cbind(as.character(DF$Species),as.character(DF$Sample))] <- DF[[mm]]
      names(output)[i] <- mm
    }
    Mat <- output
  }
  return(Mat)
}