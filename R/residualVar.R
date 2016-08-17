#' Calculates residual variance by removing 'prediction' from data matrix
#' @param prediction predicted matrix of relative abundances
#' @param Data Data matrix
#' @return Takes clr transform of Data, subtracts clr transform of prediction, and calculates the variance of all elements in the matrix.
#' @examples
#'  set.seed(1)
#' tree <- unroot(rtree(7))

#' X <- as.factor(c(rep(0,5),rep(1,5)))
#' sigClades <- Descendants(tree,c(9,12),type='tips')
#' Data <- matrix(rlnorm(70,meanlog = 8,sdlog = .5),nrow=7)
#' rownames(Data) <- tree$tip.label
#' colnames(Data) <- X
#' Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
#' Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
#' Data <- t(clo(t(Data)))

#' PF <- PhyloFactor(Data,tree,X,nclades=2)
#' prediction <- predict.phylofactor(PF,factors=1:2)
#'
#' var(c(clr(t(Data))-clr(t(prediction))))
#' residualVar(prediction,Data)

residualVar <- function(prediction,Data){
  return(sum(apply(compositions::clr(t(Data))-compositions::clr(t(prediction)),MARGIN=2,var)))
  # return(var(c(compositions::clr(t(Data))-compositions::clr(t(prediction)))))
}
