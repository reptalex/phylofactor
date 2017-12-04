#' Predicts based on phylofactorization
#'
#' @export
#' @param PF PhyloFactor object
#' @param factors Set of factors to be used for prediction. Currently, must be a continuous integer sequence from 1:nfactors
#' @param ... optional input arguments for glm
#' @return estimated data based on estimates from phylofactorization
#' @examples
#' library(caper)
#' set.seed(1)
#' tree <- unroot(rtree(8))
#' Data <- matrix(rlnorm(80),nrow=8)
#' rownames(Data) <- tree$tip.label
#' Grp <- clade.members.list(tree)$'12'
#' for (gg in Grp){
#' Data[gg,] <- Data[gg,]*20^(seq(-1,1,length.out=10))
#' }
#'
#' X <- 1:10
#' colnames(Data) <- X
#'
#' Data <- clo(Data)
#' rownames(Data) <- tree$tip.label
#'
#' PF <- PhyloFactor(Data,tree,X,nfactors=1)
#'
#' CommunityEst <- pf.predict(PF,factors=1,newdata=data.frame(X=seq(0,11,by=1/3)))
#' rownames(CommunityEst) <- rownames(Data)
#' 
#' clr <- function(A) apply(A,MARGIN=2,FUN=function(a) log(a)-mean(log(a)))
#' phytools::phylo.heatmap(tree,clr(Data))
#' phytools::phylo.heatmap(tree,clr(CommunityEst))

pf.predict <- function(PF,factors=NULL,...){
  #outputs the predictions, using nfactors from PhyloFactor.
  #additional arguments '...' are for function predict()

  if (is.null(factors)){factors=1:(PF$nfactors)}

  if (max(factors)>1 && length(factors)==1){factors=1:factors}

    coefs <- NULL
    d=1
    for (nn in factors){
      pp <- predict(PF$models[[nn]],...)
      coefs <- matrix(c(coefs,pp),ncol=d,byrow=F)
      d=d+1
    }
    
    Dat <- ilrInv(coefs,PF$basis[,factors,drop=F])
    
    rownames(Dat) <- rownames(PF$Data)
#     if (!is.null(newdata)){
#       colnames(Dat) <- colnames(PF$Data)
#     }

  return(Dat)
}
