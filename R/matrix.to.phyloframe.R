#' Converts matrix to phyloframe
#' 
#' @export
#' @param Mat Matrix. Rownames must be species and colnames must be samples found in \code{X$Sample}
#' @param X data.table or data.frame of meta-data to include in phyloframe
#' @param data.name Character. Name of data from \code{Mat} for phyloframe
#' @param empty.val Values of \code{Mat} to ignore and not include in phyloframe
#' @examples
#' library(phylofactor)
#' set.seed(1)
#' m=3   #number of species
#' n=5   #number of samples
#' d=10  #number of data points
#' DF <- data.frame('Species'=sample(letters[1:m],d,replace=T),
#'                  'Sample'=sample(1:n,d,replace=T),
#'                  'Abundance'=rpois(d,4))
#' Mat <- phyloframe.to.matrix(DF,mat.data='Abundance',empty.val=NA)
#' DF2 <- matrix.to.phyloframe(Mat,data.name='Abundance')
#' 
#' X <- data.frame('Sample'=as.character(1:n),'x'=2*(1:n),'y'=-(1:n))
#' matrix.to.phyloframe(Mat,X,data.name='Abundance')
matrix.to.phyloframe <- function(Mat,X=NULL,data.name='Data',empty.val=NA){
  if (class(Mat)!='matrix'){
    stop('input Mat must be class "matrix"')
  }
  if (!is.null(X)){
    if (!'Sample' %in% names(X)){
      stop('X must contain "Sample" to merge with phyloframe')
    }
    if (!'data.table' %in% class(X)){
      X <- data.table::as.data.table(X)
    }
    setkey(X,Sample)
  }
  species <- rownames(Mat)
  samples <- colnames(Mat)
  DF <- data.table::data.table('Species'=rep(species,times=ncol(Mat)),
                   'Sample'=rep(samples,each=nrow(Mat)),
                   'Data'=c(Mat))
  if (is.na(empty.val)){
    DF <- DF[!is.na(Data)]
  } else {
    DF <- DF[Data!=empty.val]
  }
  names(DF)[names(DF)=='Data'] <- data.name
  
  if (!is.null(X)){
    setkey(DF,Sample)
    DF <- DF[X]
  }
  return(DF)
}