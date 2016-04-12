#' Parallelized version of \code{\Link{phyloreg}} for PhyloRegression
#'
#' @param Grps Groups over which to amalgamate Data & regress
#' @param Data Data matrix whose rows are parts and columns samples
#' @param X independent variable
#' @param frmla Formula for dependence of y on x
#' @param choice Method for choosing best group. choice='var' parallelizes computation of residual variance.
#' @param method method for amalgamating
#' @param Pbasis Coming soon...
#' @param cl phyloFcluster object. See \code{\Link{phyloFcluster}}
#' @return list of glm objects
#' @example
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
parG <- lapply(clusterSplit(cl,1:m),FUN <- function(ind,g){return(g[ind])},g=Grps)

reg <- parLapply(cl,X=parG,fun=phyreg,Data=Data,XX=X,frmla=frmla,n=n,choice=choice,method=method,Pbasis=Pbasis,...)

  output <- NULL
  for (pp in 1:length(parG)){
    output$GLMs <- c(output$GLMs,reg[[pp]]$GLMs)
    output$Y <- c(output$Y,reg[[pp]]$Y)
    output$residualvar <- c(output$residualvar,reg[[pp]]$residualvar)
  }

  return(output)
}
