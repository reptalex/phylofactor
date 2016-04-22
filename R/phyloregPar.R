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

phyloregPar <- function(Grps,Data,X,frmla,choice,method,Pbasis,cl,...){

n <- dim(Data)[1]
m <- length(Grps)
parG <- lapply(parallel::clusterSplit(cl,1:m),FUN <- function(ind,g){return(g[ind])},g=Grps)
if (choice=='var'){
  reg <- parLapply(cl,parG,fun=phyreg,Data=Data,XX=X,frmla=frmla,n=n,choice=choice,method=method,Pbasis=Pbasis,...)
} else {
  # in this case, we can avoid passing the dataset to the cluster, and instead pass just the variables, Y
  Y <- Grps %>% lapply(FUN=amalgamate,Data=Data,method)
  parY <- lapply(parallel::clusterSplit(cl,1:m),FUN <- function(ind,g){return(g[ind])},g=Y)

  #X and frmlas need to be put into lists
  frmlas <- vector('list',length(parY))
  Xs <- frmlas
  for (nn in 1:length(parY)){
    frmlas[[nn]] <- frmla
    Xs[[nn]] <- X
  }
  reg <- clusterMap(cl, fun=phyreg, XX=Xs,frmla=frmlas,n=n,choice=choice,method=method,Pbasis=Pbasis, Grps=parG, Y=parY,...)

}

  output <- NULL
  for (pp in 1:length(parG)){
    # output$GLMs <- c(output$GLMs,reg[[pp]]$GLMs)
    output$stats <- rbind(output$stats,reg[[pp]]$stats)
    output$Yhat <- c(output$Yhat,reg[[pp]]$Yhat)
    output$Y <- c(output$Y,reg[[pp]]$Y)
    output$residualvar <- c(output$residualvar,reg[[pp]]$residualvar)
  }

  return(output)
}
