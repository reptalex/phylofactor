#' marginally stable aggregation of binomial data
#' 
#' @export
#' @param grp see output of \code{\link{getPhyloGroups}}
#' @param tree phylo object
#' @param Data matrix or list of matrices containing LHS variables to be aggregated
#' @param MetaData meta-data table containing "sample" and variables found in formula
#' @param expfamily Character, currently either "binomial" or "gaussian", indicating family of exponential family for marginally-stable aggregation.
#' @param frmla formula - the response variable (or cbind(Successes,Failures) for \code{expfamily='binomial'},) will be aggregated.
mAggregation <- function(Data,grp,tree,MetaData,expfamily,frmla){
  LHS <-  setdiff(as.character(frmla[[2]]),'cbind')
  if (expfamily=='binomial'){
    DF2 <- data.table::data.table('Sample'=rep(colnames(Data[[1]]),times=2),
                                  'Successes'=c(colSums(Data[['Successes']][grp[[1]],,drop=F]),colSums(Data[['Successes']][grp[[2]],,drop=F])),
                                  'Failures'=c(colSums(Data[['Failures']][grp[[1]],,drop=F]),colSums(Data[['Failures']][grp[[2]],,drop=F])),
                                  'phylo'=factor(rep(c('R','S'),each=ncol(Data[[1]]))),key='Sample')
  } else {
    DF2 <- data.table::data.table('Sample'=rep(colnames(Data),times=2),
                      'Data'=c(colSums(Data[grp[[1]],,drop=F]),colSums(Data[grp[[2]],,drop=F])),
                      'phylo'=factor(rep(c('R','S'),each=ncol(Data))),key='Sample')
    names(DF2)[2] <- LHS
  }
  DF2 <- data.table:::`[.data.table`(DF2,MetaData)
  return(DF2)
}
