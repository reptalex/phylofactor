#' Predicts based on phylofactorization
#'
#' @export
#' @param PF PhyloFactor object
#' @param factors Set of factors to be used for prediction. Currently, must be a continuous integer sequence from 1:nfactors
#' @param ... optional input arguments for glm
#' @return estimated data based on estimates from phylofactorization
#' @examples
#' library(compositions)
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
#' T <- 1:10
#' colnames(Data) <- T
#'
#' Data <- t(rcomp(t(Data)))
#' rownames(Data) <- tree$tip.label
#' phylo.heatmapAW(tree,clr(t(Data)))
#'
#' PF <- PhyloFactor(Data,tree,X=T,nclades=1)
#'
#' CommunityEst <- predict.phylofactor(PF,factors=1,newdata=data.frame(X=seq(0,11,by=1/3)))
#' rownames(CommunityEst) <- rownames(Data)
#'
#' phylo.heatmapAW(tree,clr(t(CommunityEst)))

phylofactor.predict <- function(PF,factors=NULL,...){
  #outputs the predictions, using nfactors from PhyloFactor.
  #additional arguments '...' are for function predict()

  if (is.null(factors)){factors=1: (PF$nfactors)}

    coefs <- NULL
    d=1
    for (nn in factors){
      gg <- PF$glms[[nn]]
      # pp <- predict(gg,...)
      pp <- predict(gg)
      coefs <- matrix(c(coefs,pp),ncol=d,byrow=F)
      d=d+1
    }

    Dat <- coefs %>% apply(MARGIN=1,compositions::ilrInv,V=PF$basis[,factors])

  return(Dat)
}
