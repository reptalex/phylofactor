#' objective function for generalized phylofactorization

#' @export
#' @param model model, such as glm, lm, nls, gam, which can be put into anova.
#' @param PartitioningVariables string - used to find partitioning variables in anova table.
objectiveFunction <- function(model,PartitioningVariables=''){
  
  # ss <- tryCatch(anova(model,test='Chisq'),
                 # error=function(e) 0)
  ss <- anova(model)
  if ('anova' %in% class(ss)){
    rn <- rownames(ss)
    nms <- (grepl('phylo:',rn) | grepl(':phylo',rn)) & sapply(PartitioningVariables,FUN=function(s,rn) grepl(s,rn),rn=rn)
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