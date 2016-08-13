#' Parallelized version of \code{\link{pglm}} for PhyloRegression
#'
#' @param Grps Groups over which to amalgamate Data & regress
#' @param Data Data matrix whose rows are parts and columns samples
#' @param X independent variable
#' @param frmla Formula for dependence of y on x
#' @param choice Method for choosing best group. choice='var' parallelizes computation of residual variance.
#' @param method method for amalgamating
#' @param Pbasis Coming soon...
#' @param cl phyloFcluster object. See \code{\link{phyloFcluster}}
#' @return list of glm objects
#' @examples
#' set.seed(1)
#' tree <- unroot(rtree(7))

#' X <- as.factor(c(rep(0,5),rep(1,5)))
#' sigClades <- Descendants(tree,c(9,12),type='tips')
#' Data <- matrix(rlnorm(70,meanlog = 8,sdlog = .5),nrow=7)
#' rownames(Data) <- tree$tip.label
#' colnames(Data) <- X
#' Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
#' Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
#' Data <- t(clo(t(Data)))
#' Grps <- getGroups(tree)
#' frmla <- Data ~ X
#'
#' cl <- phyloFcluster(2)
#' phyloregPar(Grps,Data,X,frmla,choice='var',method='ILR',Pbasis=1,cl=cl)
#' stopCluster(cl)
#' gc()


phyloregPar <- function(Grps,Data,X,frmla,cl,choice,Pval.Cutoff,...){

  #There are up to four tasks here that are parallelizable:
  #1) Computing ILR coordinates (amalg.ilr)
  #2) GLMs of ILR coordinates against independent variables (pglm)
  #3) calculation of Pvalues and F-statistics (getStats)
  #4) Calculation of residual variance.  (residualVar)
  # The crux here, however, is that the dataset in Data may be large, so application of parLapply, which repeatedly passes Data back and forth to clusters, will be slow.
  # For choice='F', this is not a concern: the clusters only need the full Data when calculating the residual variance.
  #Consequently, parallelization for choice='F' is done simply with parLapply in-line in PhyloRegression

  
n <- dim(Data)[1]
m <- length(Grps)
output <- NULL
output$stats <- matrix(NA,ncol=2,nrow=m)
output$Yhat <- vector(mode='list',length=m)


  parG <- lapply(parallel::clusterSplit(cl,1:m),FUN <- function(ind,g){return(g[ind])},g=Grps) #this splits our list of groups across the clusters 
  reg <- parallel::parLapply(cl,parG,fun=phyreg,Data=Data,XX=X,frmla=frmla,n=n,choice=choice,Pval.Cutoff=Pval.Cutoff,...)

  
  Ydum <- vector(mode='list',length=m)
  output$residualvar <- numeric(m)
  inds=0
  for (pp in 1:length(parG)){
    # output$GLMs <- c(output$GLMs,reg[[pp]]$GLMs)
    inds <- (max(inds)+1):(max(inds)+length(parG[[pp]]))
    output$stats[inds,] <-reg[[pp]]$stats
    output$Yhat[inds] <- reg[[pp]]$Yhat
    Ydum[inds] <- reg[[pp]]$Y
    if (choice=='var'){
    output$residualvar[inds] <- reg[[pp]]$residualvar
    }
  }
  colnames(output$stats) <- c('Pval','F')
  output$Y <- Ydum
  rm('Ydum')
  gc()
  

  return(output)
}
