#' obtains deviance for partitioning variables - default objective function for \code{\link{gpf}}

#' @export
#' @param model model, such as glm, lm, nls, gam, which can be put into anova.
#' @param grp two-member list of integers, see \code{\link{getPhyloGroups}}.
#' @param tree phylo class object
#' @param PartitioningVariables string - used to find partitioning variables in anova table.
#' @param model.fcn input model function for \code{\link{gpf}}
#' @param phyloData Data frame passed to \code{model.fcn}, used here to compute null deviance without \code{phylo} factors in gam objects
#' @param ... additional arguments passed to \code{model.fcn} passed down from the input to \code{\link{gpf}}
pvDeviance <- function(model,grp,tree,PartitioningVariables='',model.fcn,phyloData,...){
  
  if (is.null(PartitioningVariables)){
    pvs <- ''
  } else {
    pvs <- PartitioningVariables
  }
  # ss <- anova(model)
  # if (any(grepl('anova',class(ss)))){
  if ('gam' %in% class(model)){
    frmla <- model$formula
      ## EDIT: need method to extract deviance for partitioning variables in gams
      ## need likelihood of model without PVs and likelihood of models with.
      ## will have to re-run model.fcn without partitioning variables.
      ## Note: if PartitioningVariables phylo:x, we want to run a model with phylo and x, and not phylo*x
      ## likewise, for s(x,by=phylo) we want s(x)
      trms <- unique(c(attr(stats::terms(frmla),'term.labels'),setdiff(as.character(frmla[[3]]),'+')))
      
      contains.phylo <- grepl('phylo',trms)
      smoothed.by.phylo <- grepl(', by = phylo',trms)
      if (pvs==''){
        dropped.terms <- contains.phylo
      } else {
        contains.pvs <- sapply(pvs,grepl,trms)
        if (length(pvs)>1){
          contains.pvs <- apply(contains.pvs,1,any)
        }
        dropped.terms <- contains.pvs & contains.phylo
      }
      ## some of these dropped terms may be s(x,by=phylo). We don't want to drop them, but instead remove ', by=phylo'
      ix <- smoothed.by.phylo & dropped.terms
      if (any(ix)){ #a partitioning variable was smoothed by=phylo. We'll remove 'by=phylo' from these for the null model.
        trms[ix] <- gsub(', by = phylo','',trms[ix])
        ## we no longer need to drop these terms. 
        dropped.terms <- dropped.terms & !ix
      }
      trms <- paste(trms[!dropped.terms],collapse='+')
      if (model$family$family=='binomial'){
        LHS <- as.character(frmla[[2]])
        if (length(LHS)>1){ 
          if (!all.equal(c('cbind','successes','failures'),tolower(LHS))){
            stop('left-hand-side of formula for binomial must be either one entry or cbind(Successes,Failures)')
          }
          LHS <- paste(LHS[1],'(',paste(LHS[2:3],collapse=','),')',sep='')
        }
        frmla.null <- stats::as.formula(paste(LHS,'~',trms,sep=''))
      } else {
        frmla.null <- stats::as.formula(paste(frmla[[2]],'~',trms,sep=''))
      }
      
      null.model <- do.call(model.fcn,args = list('formula'=frmla.null,
                                                  'data'=phyloData[setdiff(1:nrow(phyloData),model$na.action),],...))
      
      omega <- null.model$deviance-model$deviance
      
    } else {
    
      ss <- stats::anova(model)
      if (any(grepl('anova',class(ss)))){
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
    }
  
  return(omega)
}