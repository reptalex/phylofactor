#' obtains deviance for partitioning variables - default objective function for \code{\link{gpf}}

#' @export
#' @param model model, such as glm, lm, nls, gam, which can be put into anova.
#' @param grp two-member list of integers, see \code{\link{getPhyloGroups}}.
#' @param tree phylo class object
#' @param PartitioningVariables string - used to find partitioning variables in anova table.
pvDeviance <- function(model,grp,tree,PartitioningVariables=''){
  
  if (is.null(PartitioningVariables)){
    pvs <- ''
  } else {
    pvs <- PartitioningVariables
  }
  ss <- anova(model)
  if ('anova' %in% class(ss)){
    rn <- rownames(ss)
    nms <- grepl('phylo',rn)
    rn2 <- gsub('phylo','',rn)
    vs <- sapply(pvs,FUN=function(s,rn2) grepl(s,rn2),rn2=rn2)
    if (length(pvs)>1){
      nms <- nms & apply(vs,1,any)
    } else {
      nms <- nms & vs
    }
    if (!any(nms)){
      stop('could not find objective variables and :phylo or phylo: in rownames of anova(model)')
    }
    omega <- sum(ss$Deviance[nms])
  } else {
    warning(paste('Model failed for a group. Objective set to 0 for phyloGroups of sizes',r,'and',s,sep=' '))
    omega <- 0
  }
  
  return(omega)
}