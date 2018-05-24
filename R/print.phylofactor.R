#' Prints phylofactor class objects
#' @export
#' @param PF phylofactor object
print.phylofactor <- function(PF){
  tbl <- PF$factors
  bns <- PF$bin.sizes
  grp1.sizes <- lapply(PF$groups,'[[',1) %>% sapply(length)
  biggest.grp <- which.max(grp1.sizes)
  n.singletons <- bns$Number.of.Bins[bns$Bin.Size==1]
  if (length(n.singletons)==0){
    n.singletons=0
  }
  if (PF$nfactors>10){
    tbl1 <- tbl[1:10,]
    tbl2 <- NULL
    if (biggest.grp>10){
      tbl2 <- tbl[c(biggest.grp,PF$nfactors),,drop=F]
    }
  } else {
    tbl1 <- PF$factors
    tbl2 <- NULL
  }
  if (ncol(tbl1)>2){
    tbl1[,3:ncol(tbl1)] <- apply(tbl1[,3:ncol(tbl1),drop=F],2,signif,5)
    if (!is.null(tbl2)){
      tbl2[,3:ncol(tbl2)] <- apply(tbl2[,3:ncol(tbl2),drop=F],2,signif,5)
    }
  }

  tbl.str <- paste(capture.output(print.data.frame(tbl1)),collapse='\n')
  if (!is.null(tbl2)){
  tbl.str <- paste(tbl.str,'
                      ...           ...           ...           ...           ...                              
',paste(capture.output(print.data.frame(tbl2)),collapse='\n'),sep='')
  }
  
  
  if (PF$phylofactor.fcn %in% c('PhyloFactor','gpf')){
    if (!is.null(PF$models)){
        formula <- paste('
Formula                   : ',Reduce(paste,deparse(PF$models[[1]]$formula)),sep='')
    } else {
      if (PF$choice=='custom'){
        formula <- 'Customized'
      }
      if (PF$algorithm=='CoefContrast'){
        formula <- paste('
Formula                   : ',Reduce(paste,deparse(PF$species.models[[1]]$formula)),sep='')
      }
    }
    if (PF$phylofactor.fcn=='PhyloFactor'){
      choice=paste('
Choice                    : ',PF$choice,sep='')
      algorithm=NULL
      PartitioningVariables <- NULL
      if (PF$choice=='var'){
        ExplainedVar <- paste('
% Explained Variance      : ',signif(sum(PF$factors$ExpVar),3),sep='')
      }
    } else {  # gpf
      choice=NULL
      algorithm=paste('
Algorithm                 : ',PF$algorithm,sep='')
      pvs <-PF$PartitioningVariables
      if (length(pvs)==1){
        PartitioningVariables <- paste('
Partitioning Variable     : ',pvs,sep='') 
      } else {
      pvs <- paste('{',paste(pvs,collapse=','),'}',sep='')
      PartitioningVariables <- paste('
Partitioning Variables   : ',pvs,sep='') 
      }
      ExplainedVar <- NULL
    }
  } else { #twoSamplefactor or PhyCA
    formula <- NULL
    choice <- NULL
    algorithm <- NULL
    PartitioningVariables <- NULL
    if (PF$phylofactor.fcn=='PhyCA'){
      ExplainedVar <- paste('
Frac Explained Variance   : ',signif(sum(PF$factors$ExpVar),3),sep='')
    } else {
      ExplainedVar <- NULL
    }
  }
  
  
  ln <- paste('---------------------------------',
              paste(rep('-',nchar(PF$phylofactor.fcn)),collapse=''),sep='')
  output <- paste('       phylofactor object from function ',PF$phylofactor.fcn,'
       ',ln,'       
Method                    : ',PF$method,algorithm,choice,formula,PartitioningVariables,'
Number of species         : ',length(PF$tree$tip.label),'
Number of factors         : ',PF$nfactors,ExplainedVar,'
Largest non-remainder bin : ',max(grp1.sizes),'
Number of singletons      : ',n.singletons,'
Paraphyletic Remainder    : ',length(PF$bins[[1]]),' species
                  
-------------------------------------------------------------
Factor Table:','
',tbl.str,
                  sep='')
  cat(output)
}
