#' Converts matrix to phyloframe
#' 
#' @export
#' @param Mat Matrix. Rownames must be species
#' @param MetaData data.table or data.frame of meta-data to include in phyloframe
#' @param data.name Character. Name of data from \code{Mat} for phyloframe
#' @param empty.val Values of \code{Mat} to ignore and not include in phyloframe
#' @examples
#' library(phylofactor)
#' set.seed(1)
#' m=3   #number of species
#' n=5   #number of samples
#' d=10  #number of data points
#' DF <- data.frame('Species'=sample(letters[1:m],d,replace=TRUE),
#'                  'Sample'=sample(1:n,d,replace=TRUE),
#'                  'Abundance'=rpois(d,4))
#' Mat <- phyloframe.to.matrix(DF,mat.data='Abundance',empty.val=NA)
#' DF2 <- matrix.to.phyloframe(Mat,data.name='Abundance')
#' 
#' MetaData <- data.frame('Sample'=as.character(1:n),'x'=2*(1:n),'y'=-(1:n))
#' matrix.to.phyloframe(Mat,MetaData,data.name='Abundance')
matrix.to.phyloframe <- function(Mat,MetaData=NULL,data.name='Data',empty.val=NA){
  if (class(Mat)!='matrix'){
    stop('input Mat must be class "matrix"')
  }
  if (is.null(rownames(Mat))){
    stop('Matrix must have rownames corresponding to the species')
  }
  if (!is.null(MetaData)){
    if (nrow(MetaData)!=ncol(Mat)){
      stop('number of rows in MetaData does not match the number of columns in Mat')
    }
    if (!'Sample' %in% names(MetaData)){
      if (is.null(colnames(Mat))){
        colnames(Mat) <- paste('Sample',1:ncol(Mat))
        MetaData$Sample <- paste('Sample',1:ncol(Mat))
      } else {
        MetaData$Sample <- colnames(Mat)
      }
    }
    if (!'data.table' %in% class(MetaData)){
      MetaData <- data.table::as.data.table(MetaData)
    }
    setkey(MetaData,Sample)
  }
  species <- rownames(Mat)
  samples <- colnames(Mat)
  if (is.null(samples)){
    samples <- paste('Sample',1:ncol(Mat))
  }
  DF <- data.table::data.table('Species'=rep(species,times=ncol(Mat)),
                   'Sample'=rep(samples,each=nrow(Mat)),
                   'Data'=c(Mat))
  if (is.na(empty.val)){
    DF <- DF[!is.na(Data)]
  } else {
    DF <- DF[Data!=empty.val]
  }
  names(DF)[names(DF)=='Data'] <- data.name
  
  if (!is.null(MetaData)){
    setkey(DF,Sample)
    DF <- DF[MetaData]
  }
  return(DF)
}
