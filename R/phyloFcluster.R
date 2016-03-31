phyloFcluster <- function(ncores=2,...){
  #creates cluster with all the necessary packages and functions to use these functions in parallel
  cl <- makeCluster(ncores,...)
  clusterEvalQ(cl,library(picante))
  clusterEvalQ(cl,library(ape))
  clusterEvalQ(cl,library(caper))
  clusterEvalQ(cl,library(phytools))
  clusterEvalQ(cl,library(ggtree))
  clusterEvalQ(cl,library(phangorn))
  clusterEvalQ(cl,library(compositions))
  clusterEvalQ(cl,library(stats))
  clusterExport(cl,'sourceDir')
  clusterEvalQ(cl,sourceDir('C:/Users/Big Alculus/Documents/Boulder/PhyloStats/R files/Essential PhyloFactor Files/'))
  return(cl)
}