#' objective function for \code{\link{gpf}}
#' @export
#' @param grp list containing two disjoint lists of species, such as thouse output from \code{\link{getGroups}}
#' @param tree phylo class object
#' @param DF data table containing Species, observations (N), and sample
#' @param X meta-data containing variables in formula

getObjective <- function(grp,tree,Data,X,binom.size,frmla,expfamily='gaussian',model.fcn=stats::glm,PartitioningVariables='',...){
  fit <- model.fcn(formula=frmla,data=mAggregation(Data,grp,tree,X,binom.size,expfamily),...)
  omega <- getDeviance(fit,PartitioningVariables)
  return(omega)
}