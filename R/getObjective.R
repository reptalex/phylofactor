#' objective function for \code{\link{gpf}}
#' @export
#' @param grp list containing two disjoint lists of species, such as thouse output from \code{\link{getGroups}}
#' @param tree phylo class object
#' @param Data matrix whose rows are species and columns are samples
#' @param X meta-data containing variables in formula
#' @param binom.size integer; binomial size per element of data matrix for binomial m-stable aggregation.
#' @param frmla formula for \code{model.fcn}
#' @param expfamily character string indicating manner of m-stable aggregation for \code{\link{mAggregation}}. Only "binomial" is meaningfully different.
#' @param model.fcn model function, such as \code{\link{glm}} or \code{\link{mgcv::gam}}.
#' @param PartitioningVariables character vector containing of interest for phylofactorization partitioning.
#' @param mStableAgg logical. See \code{\link{gpf}}

getObjective <- function(grp,tree,Data,X,binom.size,frmla,expfamily='gaussian',model.fcn=stats::glm,PartitioningVariables='',mStableAgg,...){
  if (mStableAgg){
    fit <- model.fcn(formula=frmla,data=mAggregation(Data,grp,tree,X,binom.size,expfamily),...)
  } else {
    fit <- model.fcn(formula=frmla,data=phyloFrame(Data,grp,tree),...)
  }
  omega <- objectiveFunction(fit,PartitioningVariables)
  return(omega)
}