#' decides Method automatically for phyca and twoSample data inputs
#' 
#' @export
#' @param Data \code{Data} input to \code{\link{PhyloFactor}}
#' @param X \code{X} input into \code{\link{PhyloFactor}}
decideMethod <- function(Data,X){
  
  if ('data.frame' %in% class(Data)){
    if (!all(c('Sample','Species','N'))){
      stop('For default decideMethod, data frame must contain c("Sample","Species","N").')
    } else {
    Data <- phylo.frame.to.matrix(DF)
    }
  }
  
  if (is.null(X)){
    if (any(1 %in% dim(Data) | is.null(dim(Data)))){
      method='twoSample'
      if (length(unique(Data))==2){
        choice='Fisher'
      } else if (ks.test(Data,rnorm(1e5,sd=sd(Data)))$p.value<0.1){
        choice='Wilcox'
      } else {
        choice='ttest'
      }
    } else {
      method='phyca'
      choice='phyca'
    }
  } else {
    stop('To avoid statistical black-box abuses, decideMethod will not choose method and choice for regression functions. Set method to contrastGLM, contrastGAM, gpf or custom and see help vignette for more details.')
  }
  return(list('method'=method,'choice'=choice))
}
