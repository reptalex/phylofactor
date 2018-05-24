#' predict newdata with phylofactor object. Only works for \code{PhyloFactor} and \code{gpf} objects.
#' @param PF phylofactor class object
#' @param newdata dataset of same form as \code{PF$Data} to be used for prediction. If \code{PF} is from \code{gpf} and \code{PF$algorithm=='mStable'}, must also input \code{newMetaData}
#' @param newMetaData required input MetaData for \code{PF$algorithm=='mStable'}.
#' @param estimate.total Logical. Whether or not to estimate column sums for \code{PhyloFactor} output. Default is false, and uses raw column sums to scale the predictions for the data.
#' @param input.total positive, numeric vector of length equal to the number of samples for \code{PhyloFactor} and \code{gpf} \code{algorithm=='mStable'} objects. Will scale predictions in each sample by the input total. 
predict.phylofactor <- function(PF,newdata=NULL,newMetaData=NULL,estimate.total=FALSE,input.total=NULL,...){
  
}