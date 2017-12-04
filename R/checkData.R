#' internal \code{PhyloFactor} quality-control function
checkData <- function(Data,tree,method){
  if ('phylo' %in% class(tree)){
    spp <- tree$tip.label
    n <- length(spp)
  } else {
    if (!all(sapply(tree,FUN=function(t) 'phylo' %in% class(t)))){
      stop('Input tree is neither phylo object nor recognizeable tree-list of phylo objects.')
    }
  }
  if (method=='twoSample'){
    if (length(Data)!=n){
      stop('Length of Data for two-sample test is not equal to number of tips in phylogeny.')
    }
    if (is.null(names(Data))){
      warning('Data for two-sample test is not named. Will assume data are in same order as tip-labels of phylogeny.')
      names(Data) <- tree$tip.label
    } else {
      Data <- Data[tree$tip.label]
    }
  } else {
    if (!all(rownames(Data) %in% tree$tip.label)){stop('some rownames of Data are not found in tree')}
    if (!all(rownames(Data)==tree$tip.label)){
      warning('rows of data are in different order of tree tip-labels - use output$data for downstream analysis, or set Data <- Data[tree$tip.label,]')
      Data <- Data[tree$tip.label,]
    }
  }
  return(Data)
}