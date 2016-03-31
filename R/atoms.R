
atoms.PhyloFactor <- function(PF,nfactors=1){
  return(atoms(PF$ilrs[,1:nfactors]))
}

################################### atoms #########################################################
atoms <- function(V=NULL,G=NULL,set=NULL){
  #input ILR sub-basis, V, or, alternatively, a list of groups G along with the set of all possible taxa, set
  # G <- Groups[Groupn-n]
  
  if (is.null(V)==F){
    G <- find.unsplit.Grps(V)
    if (is.null(dim(V))){
      set=length(V)
    } else{
      set=1:dim(V)[1]
    }
    dm <- setdiff(set,unique(unlist(G)))
    return(c(G,dm))
  } 
  
  if (is.null(G)==F && is.null(set)){
    stop('If input only G, must also include the set of all possible elements')
  } else {
    
    
    #It's possible the G has some overlap in it, so first we need to split the groups in G
    
    ## This function returns the atomic sets in overlapping vectors, x and y.
    ## Specifically, it returns the elements unique to x, unique to y, and the intersect of x,y 
    set.split <- function(x,y){
      if (length(intersect(x,y))==0){
        return(NULL)
      } else {
        C <- setdiff(union(x,y),intersect(x,y))
        return(unique(list(setdiff(x,C),setdiff(y,C),C)))
      }
    }
    
    ## the next bit of code iterates through the list of groups and, if there is overlap, splits the sets into atomic elements
    for (nn in 1:(length(G)-1)){
      for (mm in (nn+1):length(G)){
        if (any(G[[nn]] %in% G[[mm]])){
          dum <- set.split(G[[nn]],G[[mm]])
          if (length(dum)==2){
            G[c(nn,mm)] <- dum
          } else {
            G[c(nn,mm)] <- dum[1:2]
            G <- c(G,dum[3])
          }
        }
      }
    }
    
  }
  
  dm <- list(setdiff(set,unique(unlist(G))))
  return(c(G,dm))
  ### Now we need to iterate through our group and split those groups with overlap
}
