#' Produces cluster object with ncores and all necessary functions to perform phylofactorization
#' @export
#' @param ncores number of cores
#' @param ... optional input arguments for makeClsuter
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

#' frmla <- Data ~ X
#' method='ILR'
#' Grps <- getGroups(tree)
#' choice='var'
#'
#' cl <- phyloFcluster(2)
#' PhyloReg <- PhyloRegression(Data,X,frmla,Grps,method,choice,cl)
#' stopCluster(cl)
#' gc()

phyloFcluster <- function(ncores=2,...){
  #creates cluster with all the necessary packages and functions to use these functions in parallel
  cl <- parallel::makeCluster(ncores,...)
  # parallel::clusterEvalQ(cl,library(picante))
#   parallel::clusterEvalQ(cl,library(ape))
#   parallel::clusterEvalQ(cl,library(caper))
#   parallel::clusterEvalQ(cl,library(ggtree))
#   parallel::clusterEvalQ(cl,library(phangorn))
  parallel::clusterEvalQ(cl,library(compositions))
  # parallel::clusterEvalQ(cl,library(magrittr))
  parallel::clusterEvalQ(cl,library(phylofactor))
  return(cl)
}
