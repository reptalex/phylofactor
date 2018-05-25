#' predict newdata with phylofactor object. Only works for \code{PhyloFactor} and \code{gpf} objects.
#' @export
#' @param PF phylofactor class object
#' @param factor factor level to use for predictions. Will use all factors in \code{1:factor} 
#' @param newMetaData required input MetaData for \code{PF$algorithm=='mStable'}.
#' @param ... additional arguments for \code{\link{predict}}
predict.phylofactor <- function(PF,factor=NULL,frmla.phylo=NULL,newMetaData=NULL,estimate.total=FALSE,input.total=NULL,...){
  if (is.null(factor)){
    factor <- PF$nfactors
    Bins <- PF$bins
  } else {
    if (max(factor)>PF$nfactors){
      stop('factor cannot exceed PF$nfactors')
    } else {
      Bins <- bins(PF$basis[,1:factor,drop=F])
    }
  }
  if (!PF$phylofactor.fcn %in% c('PhyloFactor','gpf')){
    stop('prediction only works for regression-based phylofactorization')
  }
  
  
  
  if (PF$phylofactor.fcn=='PhyloFactor'){
    if (PF$choice=='custom'){
      stop('Cannot predict customized PhyloFactor objects')
    }
    preds <- matrix(NA,nrow=nrow(PF$Data),ncol=ncol(PF$Data))
    rownames(preds) <- rownames(PF$Data)
    colnames(preds) <- colnames(PF$Data)
    for (i in 1:length(Bins)){
      args <- c(list('data'=data.frame('Data'=colMeans(PF$transform.fcn(PF$Data[Bins[[i]],])),c(PF$X))),'formula'=PF$models[[1]]$formula,c(PF$additional.args))
      if (PF$method=='glm'){
        yhat <- do.call(stats::glm,args) %>% stats::predict(...)
      } else if (PF$method=='gam'){
        yhat <- do.call(mgcv::gam,args) %>% stats::predict(...)
      } else {
        stop('can only predict PhyloFactor objects with method="glm" or method="gam".')
      }
      preds[Bins[[i]],] <- rep(yhat,each=length(Bins[[i]]))
    }
  } else { #gpf
    
    BinMap <- data.table::data.table('Species'=PF$tree$tip.label,
                                     'phylo'=as.factor(phylobin(Bins)))
    if (PF$algorithm=='CoefContrast'){
      data.table::setkey(PF$Data,Species)
      data.table::setkey(BinMap,Species)
      Data <- data.table:::`[.data.table`(PF$Data,BinMap,nomatch=0)
      frmla <- PF$frmla
      pvs <- PF$PartitioningVariables
      RHS <- attr(stats::terms(frmla),'term.labels')
      ## there are two terms: phylo-shared coefficients in PartitioningVariables, and Species-specific coefficients.
      if (is.null(pvs) | all(sort(pvs)==sort(colnames(PF$coefficient.matrix)))){
        frmla.phylo <- update(frmla,.~phylo*.)
      } else { #pvs are not all the non-offset terms on RHS
        RHS <- unique(c(unlist(sapply(as.character(frmla[[3]]),strsplit,'\\+')),'1'))
        RHS <- gsub(' ','',RHS)
        RHS <- setdiff(RHS,'')
        non.pvs <- setdiff(RHS,pvs)
        offsets <- non.pvs[grepl('offset',non.pvs)]
        non.pvs <- setdiff(non.pvs,offsets)
        pvs <- paste(paste('phylo*',pvs,sep=''),collapse='+')
        non.pvs <- setdiff(non.pvs,'1')
        if (length(non.pvs)>0){
          non.pvs <- paste(paste('Species*',non.pvs,sep=''),collapse='+')
        }
        if ('(Intercept)' %in% pvs){
          pvs <- paste(pvs,'+phylo',sep='')
          non.pvs <- paste(non.pvs,'-Species',sep='')
        } else {
          pvs <- paste(pvs,'-phylo',sep='')
          if (length(non.pvs)>0){
            non.pvs <- paste(non.pvs,'+Species',sep='')
          } else {
            non.pvs <- 'Species'
          }
        }
        RHS <- paste(pvs,non.pvs,sep='+')
        LHS <- as.character(frmla[[2]])
        if (length(LHS)>1){
          LHS <- paste(LHS[1],'(',paste(LHS[2:3],collapse=','),')',sep='')
        }
        frmla.phylo <- stats::as.formula(paste(LHS,'~',RHS,sep=''))
        warning(paste('Did not input frmla.phylo. Will be set to',Reduce(paste,deparse(frmla.phylo))))
      }
      args <- c(list('data'=Data,'formula'=frmla.phylo),PF$additional.arguments)
      preds <- do.call(PF$model.fcn,args) %>% stats::predict(...)
    } else if (PF$algorithm=='mStable'){
      expfamily <- PF$models[[1]]$family$family
      
      
      if (expfamily=='binomial'){
        if (is.null(newMetaData)){
          preds <- matrix(NA,nrow=nrow(PF$Data[[1]]),ncol=ncol(PF$Data[[1]]))
        } else {
          preds <- matrix(NA,nrow=nrow(PF$Data[[1]]),ncol=nrow(newMetaData))
        }
      } else {
        if (is.null(newMetaData)){
          preds <- matrix(NA,nrow=nrow(PF$Data),ncol=ncol(PF$Data))  
        } else {
          preds <- matrix(NA,nrow=nrow(PF$Data),ncol=nrow(newMetaData))
        } 
      }
      rownames(preds) <- PF$tree$tip.label
      
      
      phyloData <- NULL
      for (i in 1:length(Bins)){
        grp <- list(Bins[[i]],NULL)
        pd <- mAggregation(PF$Data,grp,PF$tree,PF$MetaData,expfamily,PF$frmla.phylo)[phylo=='R']
        pd$dum <- as.integer(i)
        phyloData <- rbind(phyloData,pd)
      }
      phyloData[,phylo:=factor(dum)]
      phyloData$dum <- NULL
      
      args <- c(list('data'=phyloData,'formula'=PF$frmla.phylo),c(PF$additional.arguments))
      fit <- do.call(PF$model.fcn,args)
      
      if (!is.null(newMetaData)){
        MD <- newMetaData
        setkey(MD,'Species')
        MD <- MD[BinMap]
      } else {
        MD <- PF$MetaData
      }
      for (i in 1:length(Bins)){
        MD$phylo <- factor(rep(i,nrow(MD)))
        pp <- stats::predict(fit,newdata=MD,...)
        preds[Bins[[i]],] <- rep(pp,each=length(Bins[[i]]))
        MD$phylo <- NULL
      }
      
    } else {
      data.table::setkey(PF$Data,Species)
      data.table::setkey(BinMap,Species)
      Data <- data.table:::`[.data.table`(PF$Data,BinMap,nomatch=0)
      args <- c(list('data'=Data,'formula'=PF$frmla.phylo),c(PF$additional.arguments))
      preds <- do.call(PF$model.fcn,args) %>% predict(...)
    }
  }
  
  return(preds)
}
