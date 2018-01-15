#' Null Simulation of phylofactor object
#' 
#' @param PF class phylofactor object from \code{\link{PhyloFactor}}
#' @param reps number of replicate null simulations
#' @param nfactors number of factors to run for null simulation
#' @param nullsimFcn Either NULL, the string "shuffle", or a function taking class "phylofactor" object and producing null dataset. Default for \code{twoSampleFactor} output is to resample data with replacement. If "shuffle" will shuffle rows and columns of PF$Data. Default otherwise is to simulate standard Gaussian null data.
#' @param seed optional seed for \code{\link{set.seed}}
#' @param output output to return from each simulation. Must be in \code{names(PF)} or \code{c('ExpVar','F','All')}. If 'All', will output phylofactor object from each sim (may be memory intensive)
#' @param col.shuffle Logical. Whether or not to shuffle the columns of data. Only used if \code{nullsimFcn=='shuffle'}
#' @param row.shuffle Logical. Whether or not to shuffle the rows of data. Only used if \code{nullsimFcn=='shuffle'}
#' @param ... optional arguments to \code{\link{PhyloFactor}} or \code{\link{twoSampleFactor}}
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

pf.nullsim <- function(PF,reps,nfactors=NULL,seed=NULL,nullsimFcn=NULL,output='ExpVar',col.shuffle=T,row.shuffle=T,...){
  
  permissible.outputs <- c(names(PF),c('ExpVar','F','pvals','objective','All'))
  if (! output %in% permissible.outputs){
    stop("'output' must be in c(names(PF),c('ExpVar','F','All'))")
  }
  
  if (PF$method=='twoSample' & !output %in% c('ExpVar','pvals','objective')){
    stop('output for twoSample phylofactorizations must be either pvals or objective')
  } else {
    if (output=='ExpVar'){
      output <- 'objective'
    }
  }
  
  if (is.null(nfactors)){
    nf=PF$nfactors
  } else {
    nf=nfactors
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
  
  tm <- Sys.time()
  for (rr in 1:reps){
    if (PF$method!='twoSample'){
      if (is.null(nullsimFcn)){
        Data <- rlnorm(m*n) %>% matrix(.,nrow=m)
      } else if (nullsimFcn=='shuffle'){
        if (col.shuffle){
          if (row.shuffle){
            Data <- PF$Data[sample(1:m),sample(1:n)]
          } else {
            Data <- PF$Data[,sample(1:n)]
          }
        } else {
          if (row.shuffle){
            Data <- PF$Data[sample(1:m),]
          }
        }
      } else {
        Data <- nullsimFcn(PF)
      }
      rownames(Data) <- tree$tip.label
      
      if (rr==1){
        if (PF$method=='gpf'){
          pf <- gpf(Data,tree,PF$X,nfactors=nf,algorithm=PF$algorithm,frmla = PF$frmla,frmla.phylo=PF$frmla.phylo,PartitioningVariables = PF$PartitioningVariables,...=PF$additional.arguments)
        } else {
          pf <- PhyloFactor(Data,tree,X,nfactors=nf,method = PF$method,...)
        }
      } else {
        invisible(capture.output(pf <- PhyloFactor(Data,tree,X,nfactors=nf,method = PF$method,...)))
      }
    } else {
      if (is.null(nullsimFcn)){
        Data <- sample(PF$Data,replace=T)
        names(Data) <- tree$tip.label
      } else {
        Data <- nullsimFcn(PF)
      }
      
      if (rr==1){
        pf <- twoSampleFactor(Data,tree,nf,method = PF$choice)
      } else {
        invisible(capture.output(pf <- twoSampleFactor(Data,tree,nf,method = PF$choice)))
      }
    }
      
    if (! output == 'All'){
      if (ytype=='stat'){
        Y[[rr]] <- pf$factors[,output]
      } else {
        Y[[rr]] <- pf[output]
      }
    } else {
      Y[[rr]] <- pf
    }
    
    
    tm2 <- Sys.time()
    time.elapsed <- signif(difftime(tm2,tm,units = 'mins'),3)
    if (rr==1){
      GUI.notification <- paste('\r',rr,'null simulation completed in',time.elapsed,'minutes.   ')
    } else {
      GUI.notification <- paste('\r',rr,'null simulations completed in',time.elapsed,'minutes.    ')
    }
    
    GUI.notification <- paste(GUI.notification,'Estimated time of completion:',
                                as.character(tm+difftime(tm2,tm)*reps/rr),
                                '  \r')
    cat(GUI.notification)
    flush.console()
    
    
  }
  
  return(Y)
}