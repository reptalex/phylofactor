#' Replicate objective from phylofactor object
#' @export
#' @param grp phylo-group list
#' @param PF phylofactor object
getSignal <- function(grp,PF){
  
  if (PF$phylofactor.fcn=='PhyloFactor'){
    y <- PF$contrast.fcn(grp,PF$transform.fcn(PF$Data))
    if (PF$method=='glm'){
      model.fcn <- stats::glm
    } else if (PF$method=='gam'){
      model.fcn <- mgcv::gam
    } else if (PF$method=='max.var'){
      return(var(y))
    } else {
      model.fcn <- PF$choice.fcn
      if (!is.null(PF$choice.fcn.dependencies)){
        PF$choice.fcn.dependencies()
      }
    }
    
    if (PF$method!='custom'){
      if (PF$choice %in% c('var','F')){
        if (PF$method!='gam'){
          dataset <- PF$models[[1]]$data
          dataset$Data <- y
        } else {
          dataset <- cbind('Data'=y,PF$X)
        }
        args = c(list('data'=dataset,'formula'=PF$models[[1]]$formula),PF$additional.args)
        model <- tryCatch(do.call(model.fcn,args),error=function(e) NULL)
        if (!is.null(model)){
          omega <- getStats(model)
          if (PF$choice=='var'){
            omega <- omega['ExplainedVar']
          } else {
            omega <- omega['F']
          }
        } else {
          omega <- NA
        }
      } else {
        omega <- var(y)
      }
    } else {
      args=c(list('y'=y,'X'=PF$X,'PF.output'=FALSE),PF$additional.args)
      model <- tryCatch(do.call(model.fcn,args),error=function(e) NULL)
      if (!is.null(model)){
        omega <- model$objective
      } else {
        omega <- NA
      }
    }
  } else if (PF$phylofactor.fcn=='PhyCA'){
    y <- PF$contrast.fcn(grp,PF$transform.fcn(PF$Data))
    omega <- var(y)
  } else if (PF$phylofactor.fcn=='gpf'){
    if (PF$algorithm!='CoefContrast'){
      if (PF$algorithm=='mStable'){
        if ('family' %in% names(PF$additional.arguments)){
          family <- PF$additional.arguments$family
          if (class(family)=='function'){
            expfamily <- family()$family
          } else {
            expfamily <- family$family
          }
        } else {
          expfamily <- 'gaussian'
        }
        phyloData <- mAggregation(PF$Data,grp,PF$tree,PF$MetaData,expfamily,PF$frmla.phylo)
        args= c(list('formula'=PF$frmla.phylo,'data'=phyloData),PF$additional.arguments)
        fit <- tryCatch(do.call(PF$model.fcn,args),error=function(e) NULL)
      } else {
        phyloData <- phyloFrame(PF$Data,grp,PF$tree)
        args= c(list('formula'=PF$frmla.phylo,'data'=phyloData),PF$additional.arguments)
        fit <- tryCatch(do.call(PF$model.fcn,args),error=function(e) NULL)
      }
      if (!is.null(fit)){
        omega <- do.call(PF$objective.fcn,args=c(list('model'=fit,'grp'=grp,
                                                 'tree'=PF$tree,'PartitioningVariables'=PF$PartitioningVariables,
                                                 'model.fcn'=PF$model.fcn,'phyloData'=phyloData),PF$additional.arguments))
      } else {
        omega <- NA
      }
    } else {
      n=length(PF$tree$tip.label)
      v <- matrix(ilrvec(grp,n),ncol=n)
      if (all(PF$PartitioningVariables %in% colnames(PF$coefficient.matrix))){
        B <- PF$coefficient.matrix[,PF$PartitioningVariables,drop=F]/PF$coefficient.SE[,PF$PartitioningVariables,drop=F]
      } else {
        pvs_not_found <- PF$PartitioningVariables[!PF$PartitioningVariables %in% colnames(PF$coefficient.matrix)]
        pvs_found <- setdiff(PF$PartitioningVariables,pvs_not_found)
        ix <- sapply(pvs_not_found,FUN=function(a,b) grepl(a,b),b=colnames(PF$coefficient.matrix)) %>% apply(MARGIN=1,any)
        ix <- ix | colnames(PF$coefficient.matrix) %in% pvs_found
        if (!any(ix)){
          stop('Could neither match nor grep any PartitioningVariables in the coefficients of model. Try running stats::coefficients on your input model.fcn for a single species to determine the appropriate names for PartitioningVariables.')
        }
        B <- PF$coefficient.matrix[,ix,drop=F]/PF$coefficient.SE[,ix,drop=F]
      }
      omega <- v %*% B
      if (length(omega)>1){
        omega <- sum(omega^2)
      }
    }
    
  } else { ##twoSampleFactor
      args <- c(list('grps'=grp,'Z'=PF$Data,'p.value'=F),PF$additional.arguments)
      omega <- tryCatch(do.call(PF$model.fcn,args),error=function(e) NA)
  }
  return(omega)
}