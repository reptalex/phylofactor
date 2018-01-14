#' objective function for \code{\link{gpf}}
#' @export
#' @param grp list containing two disjoint lists of species, such as thouse output from \code{\link{getGroups}}
#' @param tree phylo class object
#' @param Data matrix whose rows are species and columns are samples
#' @param X meta-data containing variables in formula
#' @param binom.size integer; binomial size per element of data matrix for binomial m-stable aggregation.
#' @param frmla formula for \code{model.fcn}
#' @param expfamily character string indicating manner of m-stable aggregation for \code{\link{mAggregation}}. Only "binomial" is meaningfully different.
#' @param model.fcn model function, such as \code{\link{glm}} or \code{gam}.
#' @param PartitioningVariables character vector containing of interest for phylofactorization partitioning.
#' @param mStableAgg logical. See \code{\link{gpf}}
#' @param objective.fcn Objective function taking output from \code{model.fcn} as input. See \code{\link{gpf}}.
#' @param ... additional arguments for \code{model.fcn}
getObjective <- function(grp,tree,Data,X=NULL,binom.size=1,frmla,expfamily='gaussian',model.fcn=stats::glm,PartitioningVariables='',mStableAgg,objective.fcn=pvDeviance,...){
  if (mStableAgg){
    fit <- model.fcn(formula=frmla,data=mAggregation(Data,grp,tree,X,binom.size,expfamily),...)
  } else {
    fit <- model.fcn(formula=frmla,data=phyloFrame(Data,grp,tree),...)
  }
  omega <- objective.fcn(fit,grp,tree,PartitioningVariables)
  return(omega)
}