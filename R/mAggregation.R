#' marginally stable aggregation of binomial data
#' 
#' @export
#' @param grp see output of \code{\link{getGroups}}
#' @param tree phylo object
#' @param DF data table containing "Species", "sample", and "N" (N= counts of species in sample)
#' @param X meta-data table containing "sample" and variables found in formula
#' @param SampleFrame Data table containing "sample" and "G" indicating all possible combinations of samples and groups.
mAggregation <- function(DF,grp,tree,X,size=1,SampleFrame=NULL){
  if (is.null(SampleFrame)){
    SampleFrame  <- expand.grid(X$sample,c('R','S'))
    names(SampleFrame) <- c('sample','G')
    SampleFrame <- data.table::as.data.table(SampleFrame)
  }
  ff <- factorFrame(grp,tree$tip.label)
  r <- unique(ff[G=='R']$Size)
  s <- unique(ff[G=='S']$Size)
  data.table::setkey(DF,Species)
  DF2 <- DF[ff]
  # DF2 <- data.table:::`[.data.table`(DF,ff)
  setkey(DF2,G,sample)
  
  ###################### marginal-stable aggregatoin ######
  ## insert code to add zeros here:
  DF2 <- DF2[,list(Successes=sum(N),
                   Failures=size*(unique(Size))-sum(N),
                   Size=unique(Size)),by=list(G,sample)]
  setkey(DF2,sample)
  
  ##### Inputting zeros
  setkey(SampleFrame,sample,G)
  setkey(DF2,sample,G)
  DF2 <- DF2[SampleFrame]
  ix <- is.na(DF2$Successes)
  DF2[ix & G=='R']$Size <- r
  DF2[ix & G=='S']$Size <- s
  DF2[ix & G=='R']$Failures <- size*r
  DF2[ix & G=='S']$Failures <- size*s
  DF2[ix]$Successes <- 0
  
  setkey(X,sample)
  DF2 <- DF2[X]
  return(DF2)
}