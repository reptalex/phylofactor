#' Get list of groups defined by phylogeny
#' @export
#' @param tree Phylogeny
#' @return list of groups and their complements, i.e. lists of groups split by each edge of the unrooted phylogeny
#' @examples
#' tr <- rtree(8)
#' tr$tip.label <- 1:8
#' par(mfrow=c(1,1))
#' plot.phylo(tr,use.edge.length=F)
#' Grps <- getGroups(tr)

######################### getGroups ######################################

getGroups <- function(tree){
  set=1:length(tree$tip.label)
  n=length(set)
  Grps <- vector(mode='list',length=(2*n-3))

  Grps[1:n] <- lapply(as.list(set[1:n]),FUN = function(x,set){return(list(x,setdiff(set,x)))},set=set)
  names(Grps)[1:n] <- 1:n

  if (n>3){
    cml <- caper::clade.members.list(tree)
    cml <- cml[2:length(cml)]
    cml <- cml[which(sapply(cml,FUN= function(x,n) length(x)<(n-1),n=n,simplify=T))]
    if (length(cml)>1){
      Grps[(n+1):(n+length(cml))] <- lapply(cml,FUN = function(x,set){return(list(x,setdiff(set,x)))},set=set)
      names(Grps)[(n+1):(n+length(cml))] <- names(cml)
    }
  }
  lapply(Grps,FUN=function(x) lapply(x,sort))
  return(Grps)
}
