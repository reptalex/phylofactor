#' Predicts based on phylofactorization
#'
#' @export
#' @param PF PhyloFactor object
#' @param factors Set of factors to be used for prediction. Currently, must be a continuous integer sequence from 1:nfactors
#' @param ... optional input arguments for predict
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
  if (PF$method!='gpf'){
    if (max(factors)>1 && length(factors)==1){factors=1:factors}
  
      coefs <- NULL
      d=1
      for (nn in factors){
        pp <- stats::predict(PF$models[[nn]],...)
        coefs <- matrix(c(coefs,pp),ncol=d,byrow=F)
        d=d+1
      }
      
      Dat <- ilrInv(coefs,PF$basis[,factors,drop=F])
      
      rownames(Dat) <- rownames(PF$Data)
  } else {
    frmla <- PF$models[[1]]$formula
    family <- PF$models[[1]]$family
    if (family$family!='binomial'){
      LHS <- as.character(frmla[[2]])
      if (pf$algorithm=='mStable'){
        m=nrow(PF$Data)
        n=ncol(PF$Data)
        samples <- colnames(PF$Data)
        Dat <- data.table('Data'=c(PF$Data),
                          'Sample'=rep(samples,each=m),
                          'phylo'=rep(factor(phylobin(PF$bins)),times=n),
                          key='Sample')[PF$MetaData] 
        names(Dat)[1] <- LHS
        Dat <- PF$model.fcn(frmla,family=family,data=Dat) %>%
          stats::predict(.,...) %>%
          matrix(.,nrow=m,ncol=n,byrow=F)
        colnames(Dat) <- colnames(PF$Data)
        rownames(Dat) <- rownames(PF$Data)
      } else {
        pbin <- data.frame('Species'=PF$tree$tip.label,
                           'phylo'=factor(phylobin(PF$bins)))
        setkey(PF$Data,'Species')
        Dat <- PF$Data[pbin] %>%
                  PF$model.fcn(frmla,family=family,data=.) %>%
                  stats::predict(.,...)
      }
    } else {
      if (pf$algorithm=='mStable'){
        m=nrow(PF$Data[[1]])
        n=ncol(PF$Data[[1]])
        samples <- colnames(PF$Data[[1]])
        Dat <- data.table('Successes'=c(PF$Data[[1]]),
                        'Failures'=c(PF$Data[[2]]),
                        'Sample'=rep(samples,each=m),
                        'phylo'=rep(factor(phylobin(PF$bins)),times=n),
                        key='Sample')[PF$MetaData] %>%
        PF$model.fcn(frmla,family=family,data=.) %>%
        stats::predict(.,...) %>%
        matrix(.,nrow=m,ncol=n,byrow=F)
        colnames(Dat) <- colnames(PF$Data[[1]])
        rownames(Dat) <- rownames(PF$Data[[1]])
      } else {
        pbin <- data.frame('Species'=PF$tree$tip.label,
                           'phylo'=factor(phylobin(PF$bins)))
        setkey(PF$Data,'Species')
        Dat <- PF$Data[pbin] %>%
          PF$model.fcn(frmla,family=family,data=.) %>%
          stats::predict(.,...)
      }
    }

                     
  }
  return(Dat)
}
