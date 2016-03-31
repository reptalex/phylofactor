######################### getGroups ######################################

getGroups <- function(tree){
  set=1:length(tree$tip.label)
  n=length(set)
  Grps <- vector(mode='list',length=(2*n-3))
  
  Grps[1:n] <- lapply(as.list(set[1:n]),FUN = function(x,set){return(list(x,setdiff(set,x)))},set=set)
  names(Grps)[1:n] <- 1:n
  
  
  cml <- clade.members.list(tree)[2:(n-2)]
  Grps[(n+1):(2*n-3)] <- lapply(cml,FUN = function(x,set){return(list(x,setdiff(set,x)))},set=set)
  names(Grps)[(n+1):(2*n-3)] <- names(cml)
  
  return(Grps)
}
