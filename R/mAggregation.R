#' marginally stable aggregation of binomial data
#' 
#' @export
#' @param grp see output of \code{\link{getGroups}}
#' @param tree phylo object
#' @param DF data table containing "Species", "sample", and "N" (N= counts of species in sample)
#' @param X meta-data table containing "sample" and variables found in formula
#' @param binom.size Size of binomial (binom.size=1 for Bernoulli)
#' @param expfamily Character, currently either "binomial" or "gaussian", indicating family of exponential family for marginally-stable aggregation.
mAggregation <- function(DF,grp,tree,X,binom.size,expfamily){
  if (expfamily=='binomial'){
    r <- length(grp[[1]])
    s <- length(grp[[2]])
    DF2 <- data.table('Sample'=rep(colnames(DF),times=2),
                      'Successes'=c(colSums(DF[grp[[1]],,drop=F]),colSums(DF[grp[[2]],,drop=F])),
                      'Failures'=rep(c(binom.size*r,binom.size*s),each=ncol(DF)),
                      'phylo'=factor(rep(c('R','S'),each=ncol(DF))),key='Sample')
  } else {
    DF2 <- data.table('Sample'=rep(colnames(DF),times=2),
                      'Data'=c(colSums(DF[grp[[1]],,drop=F]),colSums(DF[grp[[2]],,drop=F])),
                      'phylo'=factor(rep(c('R','S'),each=ncol(DF))),key='Sample')
  }
  DF2 <- DF2[X]
  return(DF2)
}