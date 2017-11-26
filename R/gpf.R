#' Generalized phylofactorization - currently skeleton for binomial glm
#' 
#' @export
#' @param DF data table containing columns of "Species", 'N' (counts, 0<=N<=3), "sample". 
#' @param X meta-data containing "sample" and variables found in \code{frmla}
#' @param nfactors integer for number of factors to find
#' @param ncores integers for number of cores to use for parallelization
#' @param size size of binomial samples for each species. Size=1 for presence/absence data.
gpf <- function(DF,X,tree,frmla,nfactors,ncores=NULL,size=1){
  if (!any(class(DF)!='data.table')){
    DF <- data.table::as.data.table(DF)
  }
  if (!any(class(X)!='data.table')){
    X <- data.table::as.data.table(X)
  }
  if (!all(c('Species','sample') %in% names(DF))){
    stop('DF must contain columns with "Species" and "sample"')
  } else {
    if (!all(unique(DF$Species) %in% tree$tip.label)){
      stop('Species in DF not found in tree')
    }
  }
  
  SampleFrame  <- expand.grid(X$sample,c('R','S'))
  names(SampleFrame) <- c('sample','G')
  SampleFrame <- data.table::as.data.table(SampleFrame)
  
  
  if (!is.null(ncores)){
    cl <- phyloFcluster(ncores)
    parallel::clusterEvalQ(cl,library(data.table))
    parallel::clusterExport(cl,varlist = c('DF','X','tree','SampleFrame'),envir = environment())
  } else {
    cl <- NULL
  }
  
  
  
  
  
  treeList <- list(tree)
  binList <- list(1:ape::Ntip(tree))
  Grps=getGroups(tree)
  output <- NULL
  pfs=0
  while (pfs < nfactors){
    
    if (pfs>=1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      Grps <- getNewGroups(tree,treeList,binList)
    }
    
    if (is.null(ncores)){
      obj <- sapply(Grps,getObjective,tree,DF,X,SampleFrame,size,frmla)
    } else {
      obj <- parallel::parSapply(cl,Grps,FUN=function(grp,tree,DF,xx,SampleFrame,size,frmla) getObjective(grp,tree,DF,xx,SampleFrame,size,frmla),tree=tree,DF=DF,xx=X,SampleFrame=SampleFrame,size=size,frmla=frmla)
    }
    
    winner <- which.max(obj)
    grp <- getLabelledGrp(tree=tree,Groups=Grps[[winner]])
    output$groups <- c(output$groups,list(Grps[[winner]]))
    
    grpInfo <- matrix(c(names(grp)),nrow=2)
    output$factors <- cbind(output$factors,grpInfo)
    
    pfs=pfs+1
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
  output$Data <- DF
  output$method <- 'gpf'
  class(output) <- 'phylofactor'
  
  return(output)
}