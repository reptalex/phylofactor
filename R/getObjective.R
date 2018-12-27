#' objective function for \code{\link{gpf}}
#' @export
#' @param grp list containing two disjoint lists of species, such as thouse output from \code{\link{getGroups}}
#' @param tree phylo class object
#' @param Data If \code{mStableAgg==TRUE}, a matrix whose rows are species and columns are samples. Otherwise, a data table whose columns include "Species" and "Sample" and whose key is "Species".
#' @param frmla formula for \code{model.fcn}
#' @param MetaData meta-data containing variables in formula and the column "Sample". If \code{mStableAgg==F}, this input is not used - all variables must be contained in \code{Data}
#' @param PartitioningVariables character vector containing of interest for phylofactorization partitioning.
#' @param mStableAgg logical. See \code{\link{gpf}}
#' @param expfamily character string indicating manner of m-stable aggregation for \code{\link{mAggregation}}. Only "binomial" is meaningfully different.
#' @param model.fcn model function, such as \code{\link{glm}} or \code{gam}.
#' @param objective.fcn Objective function taking output from \code{model.fcn} as input. See \code{\link{gpf}}.
#' @param ... additional arguments for \code{model.fcn}
getObjective <- function(grp,tree,Data,frmla,MetaData=NULL,PartitioningVariables='',mStableAgg,expfamily='gaussian',model.fcn=stats::glm,objective.fcn=pvDeviance,ignore.tips=F,...){
  if (ignore.tips){
    
    if (any(sapply(grp,length)==1)){
      omega <- -Inf
    } else {
      if (mStableAgg==T){
        phyloData <- mAggregation(Data,grp,tree,MetaData,expfamily,frmla)
        fit <- do.call(model.fcn,args=list('formula'=frmla,'data'=phyloData,...))
      } else {
        phyloData <- phyloFrame(Data,grp,tree)
        fit <- do.call(model.fcn,args=list('formula'=frmla,'data'=phyloData,...))
      }
      
      omega <- objective.fcn(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...)
    }
  } else {
  
    if (mStableAgg==T){
      phyloData <- mAggregation(Data,grp,tree,MetaData,expfamily,frmla)
      fit <- do.call(model.fcn,args=list('formula'=frmla,'data'=phyloData,...))
    } else {
      phyloData <- phyloFrame(Data,grp,tree)
      fit <- do.call(model.fcn,args=list('formula'=frmla,'data'=phyloData,...))
    }
    
    omega <- objective.fcn(fit,grp,tree,PartitioningVariables,model.fcn,phyloData,...)
  }
  return(omega)
}
