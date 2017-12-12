#' Generalized phylofactorization - currently skeleton for binomial glm; will expand to exponential, gamma etc. soon.
#' 
#' @export
#' @param Data data table containing columns of "Species", 'N' (counts, 0<=N<=3), "Sample". 
#' @param X meta-data containing "Sample" and variables found in \code{frmla}
#' @param frmla Formula. If \code{expfamily='binomial'}, must have c(Successes,Failures)~. Otherwise, the variable for data is "Data", e.g. \code{Data~phylo}. Explanatory variables used for phylofactorization must interact with the variable \code{phylo}, e.g. partitioning by pH and not N requires \code{Data~N+phylo*pH}
#' @param nfactors integer for number of factors to find
#' @param ncores integers for number of cores to use for parallelization
#' @param binom.size binom.size of binomial samples for each species. binom.size=1 for presence/absence data.
#' @param expfamily Either "gaussian" or "binomial" - determines the aggregation method.
#' @param model.fcn Regression function, such as glm, gam, glm.nb, gls. Must have column labelled "Deviance" in \code{\link{anova}}.
#' @param PartitioningVariables Character vector containing the variables in \code{frmla} to be used for phylofactorization. Objective function will be the sum of deviance from all variables listed here.
#' @examples 
#' library(phylofactor)
#' 
#' set.seed(1)
#' m <- 50
#' n <- 200
#' tree <- rtree(m)
#' X <- data.table('y'=rnorm(n),
#'                 'z'=rnorm(n,sd=0.5),
#'                 'Sample'=sapply(1:n,FUN=function(s) paste('Sample',s)),
#'                 key='Sample')
#' #we'll partition by 'y'.
#' binom.size=3
#' clade <- phangorn::Descendants(tree,75,'tips')[[1]]
#' clade2 <- phangorn::Descendants(tree,53,'tips')[[1]]
#' 
#' ######## presence/absence dataset with affected clade #######
#' ## most species have higher P{present} with y
#' ilogit <- function(eta) 1/(1+exp(-eta))
#' eta <- X$z+X$y
#' p <- ilogit(eta)
#' M <- matrix(rbinom(m*n,binom.size,rep(p,times=m)),nrow=m,ncol=n,byrow=T)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- X$Sample
#' 
#' #### the first clade decreases with y ####
#' eta1 <- X$z-X$y
#' p1 <- ilogit(eta1)
#' for (species in clade){
#'    M[species,] <- rbinom(n,binom.size,p1)
#' }
#' #### the second clade weakly decreases with y ####
#' eta2 <- X$z-.3*X$y
#' p2 <- ilogit(eta2)
#' for (species in clade2){
#'    M[species,] <- rbinom(n,binom.size,p2)
#' }
#' 
#' ####################### to partition on y, must have phylo* #########
#' pf <- gpf(M,tree,X,frmla=cbind(Successes,Failures)~z+phylo*y,nfactors=2,binom.size=binom.size,family=binomial(link='logit'))
#' all.equal(pf$groups[[1]][[1]],clade)
#' 
#' pf.tree(pf)
#' par(mfrow=c(2,1))
#' phytools::phylo.heatmap(tree,M[,order(X$y)])
#' phytools::phylo.heatmap(tree,ilogit(pf.predict(pf)[,order(X$y)]))
#' 
#' ################# Poisson Regression
#' eta <- .3*X$z
#' lambda <- exp(eta)
#' M <- matrix(rpois(m*n,rep(lambda,times=m)),nrow=m,ncol=n,byrow=T)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- X$Sample
#' 
#' #### the first clade decreases with y ####
#' eta1 <- .3*X$z-X$y
#' lambda1 <- exp(eta1)
#' for (species in clade){
#'    M[species,] <- rpois(n,lambda1)
#' }
#' #### the second clade strongly increases with y ####
#' eta2 <- .3*X$z-0.3*X$y
#' lambda2 <- exp(eta2)
#' for (species in clade2){
#'    M[species,] <- rpois(n,lambda2)
#' }
#' 
#' 
#' 
#' ##For non-binomial, use "Data" as response variable #########
#' pf <- gpf(M,tree,X,frmla=Data~z+phylo*y,nfactors=2,family=poisson)
#' list(clade,clade2)
#' pf$bins
#' 
#' par(mfrow=c(2,1))
#' phytools::phylo.heatmap(tree,M[,order(X$y)])
#' phytools::phylo.heatmap(tree,exp(pf.predict(pf)[,order(X$y)]))
gpf <- function(Data,tree,X,frmla,nfactors=NULL,ncores=NULL,binom.size=1,expfamily='gaussian',model.fcn=stats::glm,PartitioningVariables='',...){
  
  output <- NULL
  output$call <-match.call()
  output$Data <- Data
  output$X <- X
  output$model.fcn <- model.fcn
  output$additional.arguments <- list(...)
  output$binom.size=binom.size
  
  ############## Checking Data and X compatibility ####################
  if (!'data.table'%in%class(X)){
    X <- data.table::as.data.table(X)
    if (!'Sample'%in%names(X)){
      stop('meta-data X must contain a column named "Sample"')
    }
    setkey(X,Sample)
  }
  if ('data.table' %in% class(Data)){
    if (!all(c('Species','Sample') %in% names(Data))){
      stop('Data must contain columns with "Species" and "Sample"')
    } else {
      if (!all(unique(Data$Species) %in% tree$tip.label)){
        stop('Species in Data not found in tree')
      }
    }
    if (!any(X$Sample %in% Data$Sample)){
      stop('No X$Sample in Data$Sample')
    } else if (!all(X$Sample %in% Data$Sample)){
      warning(paste('Only',length(intersect(X$Sample,Data$Sample)),'samples shared in both Data and X'))
    }
    Data <- phylo.frame.to.matrix(DF)
  }
  if ('matrix' %in% class(Data)){
    if (!all(tree$tip.label==rownames(Data))){
      if (!all(tree$tip.label %in% rownames(Data))){
        stop('Not all tree tip labels are in data')
      } else if (!all(rownames(Data) %in% tree$tip.label)){
        stop('Not all rownames of data are in tree tip labels')
      } else {
        Data <- Data[tree$tip.label,]
      }
    }
    if (!any(X$Sample %in% colnames(Data))){
      stop('No X$Sample in colnames(Data)')
    } else if (!all(X$Sample %in% colnames(Data))){
      warning(paste('Only',length(intersect(X$Sample,colnames(Data))),'samples shared in both Data and X'))
    }
  }
  
  if ('family' %in% names(list(...))){
    family <- list(...)$family
    if (class(family)=='function'){
      expfamily <- family()$family
    } else {
      expfamily <- family$family
    }
  }
  
  if (!is.null(ncores)){
    cl <- phyloFcluster(ncores)
    parallel::clusterEvalQ(cl,library(data.table))
    parallel::clusterExport(cl,varlist = c('Data','X','tree','binom.size','model.fcn'),envir = environment())
  } else {
    cl <- NULL
  }
  
  
  
  
  
  treeList <- list(tree)
  binList <- list(1:ape::Ntip(tree))
  Grps=getPhyloGroups(tree)
  pfs=0
  tm <- Sys.time()
  while (pfs<nfactors){
    
    if (pfs>=1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      Grps <- getNewGroups(tree,treeList,binList)
    }
    
    if (is.null(ncores)){
      obj <- sapply(Grps,getObjective,tree,Data,X,binom.size,frmla,expfamily,model.fcn,PartitioningVariables,...)
      # obj <- sapply(Grps,getObjective,tree,Data,X,binom.size,frmla,expfamily,model.fcn,PartitioningVariables,family=family)
      
    } else {
      obj <- parallel::parSapply(cl,Grps,FUN=function(grp,tree,Data,xx,binom.size,frmla,expfamily,model.fcn,PartitioningVariables,...) getObjective(grp,tree,Data,xx,binom.size,frmla,expfamily,model.fcn,PartitioningVariables,...),tree=tree,Data=Data,xx=X,binom.size=binom.size,frmla=frmla,expfamily=expfamily,model.fcn=model.fcn,PartitioningVariables=PartitioningVariables,...)
      # obj <- parallel::parSapply(cl,Grps,FUN=function(grp,tree,Data,xx,binom.size,frmla,expfamily,model.fcn,PartitioningVariables) getObjective(grp,tree,Data,xx,binom.size,frmla,expfamily,model.fcn,PartitioningVariables),tree=tree,Data=Data,xx=X,binom.size=binom.size,frmla=frmla,expfamily=expfamily,model.fcn=model.fcn,PartitioningVariables=PartitioningVariables)
      
    }
    
    winner <- which.max(obj)
    grp <- Grps[[winner]]
    if (pfs==0){
      output$models <- list(model.fcn(frmla,data=mAggregation(Data,grp,tree,X,binom.size,expfamily),...))
    } else {
      output$models <- c(output$models,list(model.fcn(frmla,data=mAggregation(Data,grp,tree,X,binom.size,expfamily),...)))
    }
    grp <- getLabelledGrp(tree=tree,Groups=Grps[[winner]])
    output$groups <- c(output$groups,list(Grps[[winner]]))
    
    grpInfo <- matrix(c(names(grp)),nrow=2)
    output$factors <- cbind(output$factors,grpInfo)
    
    pfs=pfs+1
    tm2 <- Sys.time()
    time.elapsed <- signif(difftime(tm2,tm,units = 'mins'),3)
    if (pfs==1){
      GUI.notification <- paste('\r',pfs,'factor completed in',time.elapsed,'minutes.   ')
    } else {
      GUI.notification <- paste('\r',pfs,'factors completed in',time.elapsed,'minutes.    ')
    }
    if (!is.null(nfactors)){
      GUI.notification <- paste(GUI.notification,'Estimated time of completion:',
                                as.character(tm+difftime(tm2,tm)*nfactors/pfs),
                                '  \r')
    }
    cat(GUI.notification)
    flush.console()
  }
  
  if (!is.null(ncores)){
    parallel::stopCluster(cl)
    rm('cl')
  }
  
  output$factors <- t(output$factors)
  
  output$basis <- matrix(NA,nrow=ape::Ntip(tree),ncol=pfs)
  for (i in 1:length(output$groups)){
    output$basis[,i] <- ilrvec(output$groups[[i]],ape::Ntip(tree))
  }
  
  output$bins <- bins(output$basis)
  output$tree <- tree
  output$nfactors <- pfs
  output$Data <- Data
  output$method <- 'gpf'
  class(output) <- 'phylofactor'
  
  return(output)
}