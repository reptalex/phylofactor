#' Null Simulation of phylofactor object
#' 
#' @export
#' @param PF class phylofactor object from \code{\link{PhyloFactor}}
#' @param reps number of replicate null simulations
#' @param seed optional seed for \code{\link{set.seed}}
#' @param method string, either "Gaussian" or "Shuffle". If "Gaussian", simulations will be standard log-normal datasets. Otherwise, simulations will be original data with reshuffled rows and columns
#' @param output output to return from each simulation. Must be in \code{names(PF)} or \code{c('ExpVar','F','All')}. If 'All', will output phylofactor object from each sim (may be memory intensive)
#' @param ... optional arguments to \code{\link{PhyloFactor}}
#' @examples 
#' library(ape)
#' library(phylofactor)
#' set.seed(1)
#' m=7
#' n=10
#' tree <- rtree(m)
#' X <- rnorm(n)
#' Data <- rlnorm(m*n) %>% matrix(.,nrow=m)
#' rownames(Data) <- tree$tip.label
#' clade <- c('t1','t3','t6')
#' for (tip in clade){ Data[tip,] <- Data[tip,]*exp(8*X) }
#' PF <- PhyloFactor(Data,tree,X,nfactors=4)
#' 
#' nullsim <- pf.nullsim(PF,10,nfactors=4)
#' 
#' plot(PF$ExplainedVar,type='l')
#' for (nn in 1:10){lines(nullsim[[nn]],col='grey')}
#' legend('center',c('Original Data','Null Data'),col=c('black','gray'),lty=c(1,1))

pf.nullsim <- function(PF,reps,seed=NULL,method='Gaussian',output='ExpVar',...){
  if (! method %in% c('Gaussian','Shuffle')){
    stop('input "method" must be either "Gaussian" or "Shuffle"')
  }
  
  permissible.outputs <- c(names(PF),c('ExpVar','F','All'))
  if (! output %in% permissible.outputs){
    stop("'output' must be in c(names(PF),c('ExpVar','F','All'))")
  }
  
  if (output %in% c('ExpVar','F')){
    ytype <- 'stat'
  } else {
    ytype <- 'list'
  }
  Y <- vector(mode='list',length=reps)
  
  if (!is.null(seed)){
    base::set.seed(seed)
  }
  
  Data <- PF$Data
  tree <- PF$tree
  X <- PF$X
  m <- nrow(Data)
  n <- ncol(Data)
  
  for (rr in 1:reps){
    if (method=='Gaussian'){
      Data <- rlnorm(m*n) %>% matrix(.,nrow=m)
    } else {
      Data <- PF$Data[sample(1:m),sample(1:n)]
    }
    rownames(Data) <- tree$tip.label
    
    pf <- PhyloFactor(Data,tree,X,nfactors=pf$nfactors,...)
      
    if (! output == 'All'){
      if (ytype=='stat'){
        Y[[rr]] <- pf$factors[,output]
      } else {
        Y[[rr]] <- pf[output]
      }
    } else {
      Y[[rr]] <- pf
    }
  }
  
  return(Y)
}