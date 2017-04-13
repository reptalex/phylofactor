#' Simulates brownian motion with HGT between random clades via poisson pt process
#' 
#' @export
#' @param tree phylo class object
#' @param rate positive real number indicating rate of horizontal gene transfer per unit edge-length in tree
#' @param a real number ancestral state
#' @param mu drift
#' @param sig2 volatility
#' @examples 
#' library(phylofactor)
#' library(ape)
#' 
#' tree <- rtree(10)
#' sim <- simHGT(tree,0.01,2,-3,0.5)

simHGT <- function(tree,rate,a=0,mu=0,sig2=1){
  # Gillespie simulation of 
  n <- length(tree$tip.label)
  tipheights <- diag(ape::vcv.phylo(tree))
  Tmax <- max(tipheights)
  
  nodes <- ape::Ntip(tree)+2:ape::Nnode(tree)
  grps <- sapply(as.list(nodes),FUN=function(n,t) phangorn::Descendants(t,n,type='tips'),t=tree)
  ndheights <- sapply(as.list(nodes),FUN=function(n,t) nodeheight(t,n),t=tree)
  names(ndheights) <- nodes
  names(grps) <- nodes
  
  approx.hgt <- ceiling(rate*sum(tree$edge)) #approximate number of HGT events we will see
  
  X <- matrix(NA,ncol=(2*n-1)+2*approx.hgt,nrow=n)
  rownames(X) <- tree$tip.label
  X[,1] <- a
  binder <- matrix(NA,ncol=1,nrow=nrow(X)) ## This is used to grow X as needed
  
  tt <- 0     #stopwatch
  event <- 1  #event counter
  hgt <- 0   #hgt event counter
  tvec <- 0   #vector of event times, will be equal to ncol(X).
  HGTs <- matrix(NA,nrow=3,ncol=2*approx.hgt)
  nodes.in.play <- phangorn::Descendants(tree,n+1,type='children') 
  nn <- length(nodes.in.play)  # Number of edges currently evolving
  nt <- length(tipheights)     # Number of tips still in play
  
  while (tt<Tmax){
    
    ## The next event is either: HGT, a node, or a tip.
    ## First, find distance to next tip or node
    possibilities <- c(tipheights,ndheights)
    ix.min <- which.min(possibilities)
    dtmax <- possibilities[ix.min]-tt
    n.hgts <- rpois(1,lambda=rate*dtmax*nn) #Number of HGT events between now and next split/tip

    ######### INCOMPLETE - calculate states (X+dx), HGTs, (X[ix1]=X[ix2]) and times (t+dt)
    if (n.hgts>0){
      
      ts <- runif(n.hgts,min=tt,max=tt+dtmax) %>% sort
      for (i in 1:n.hgts){
        event=event+1
        hgt=hgt+1
        
        
        ####### Update time
        if (i==1){
          DT <- ts[i]-tt
        } else {
          DT <- ts[i]-ts[i-1]
        }
        tvec <- c(tvec,tvec[event-1]+DT)
        
        ####### Update X
        if (event>ncol(X)){
          X <- cbind(X,binder)
        }
        
        dx <- rnorm(nn,mean=mu*DT,sd = sqrt(sig2*DT)) #change happens before HGT
        for (ll in 1:nn){
          nde <- nodes.in.play[ll]
          if (nde>n){
            ixs <- tree$tip.label[grps[toString(nde)][[1]]]
          } else {
            ixs <- tree$tip.label[nde]
          }
          X[ixs,event] <- X[ixs,event-1]+dx[ll]
        }
        
        
        ### Perform HGT. Randomly sample receiver and donor node
        if (hgt>ncol(HGTs)){
          hb <- c(sample(nodes.in.play,2,replace=F),tvec[event])
          HGTs <- cbind(HGTs,matrix(hb,ncol=1))
        } else {
          HGTs[1:2,hgt] <- sample(nodes.in.play,2,replace=F)
          HGTs[3,hgt] <- tvec[event]
        }
        
        
        if (HGTs[1,hgt]<=n){ # receiving node is a tip
          ixr <- tree$tip.label[HGTs[1,hgt]]
        } else {
          ixr <- tree$tip.label[grps[toString(HGTs[1,hgt])][[1]]] #receiving node
        }
        if (HGTs[2,hgt]<=n){ # receiving node is a tip
          ixd <- tree$tip.label[HGTs[2,hgt]]
        } else {
          ixd <- tree$tip.label[grps[toString(HGTs[2,hgt])][[1]]] #receiving node
        } #donor node
        
        X[ixr,event] <- X[ixd[1],event]  ### replace receiver traits with donor traits
        ixn <- setdiff(rownames(X),c(ixr,ixd))
        X[ixn,event] <- X[ixn,event-1]
      }
      
      
      ## simulate the last bit of drift until our next split or terminal node
      event=event+1
      if (event>ncol(X)){
        X <- cbind(X,binder)
      }
      tvec <- c(tvec,tt+dtmax)
      DT <- tvec[length(tvec)]-tvec[length(tvec)-1]
      dx <- rnorm(nn,mean=mu*DT,sd = sqrt(sig2*DT)) #change happens before HGT
      for (ll in 1:nn){
        nde <- nodes.in.play[ll]
        if (nde>n){
          ixs <- tree$tip.label[grps[toString(nde)][[1]]]
        } else {
          ixs <- tree$tip.label[nde]
        }
        X[ixs,event] <- X[ixs,event-1]+dx[ll]
      }
      
    } else { ### no HGT
      
      event=event+1
      if (event>ncol(X)){
        X <- cbind(X,binder)
      }
      tvec <- c(tvec,tt+dtmax)
      dx <- rnorm(nn,mean=mu*dtmax,sd = sqrt(sig2*dtmax))
      
      for (ll in 1:length(nodes.in.play)){
        nde <- nodes.in.play[ll]
        if (nde>n){
          ixs <- tree$tip.label[grps[toString(nde)][[1]]]
        } else {
          ixs <- tree$tip.label[nde]
        }
        X[ixs,event] <- X[ixs,event-1]+dx[ll]
      }
      
    }
    
    
    ######### INCOMPLETE - Update. Remove tips, add nodes, change nodes.in.play
    if (ix.min <= nt){ #we finish at a tip. Need to remove it from the pool
      tipheights <- setdiff(tipheights,tipheights[ix.min])
      nt=nt-1
    } else { #we finish at a node - need to remove it from the pool and get its decsendants
      rm.node <- names(ndheights)[ix.min-nt]
      add.nodes <- phangorn::Descendants(tree,rm.node,type='children')
      nodes.in.play <- c(setdiff(nodes.in.play,rm.node),add.nodes)
      tips <- tree$tip.label[add.nodes[which(add.nodes<=n)]]
      add.nodes <- add.nodes[which(add.nodes>n)]
      ndheights <- ndheights[setdiff(names(ndheights),rm.node)]
      nn <- length(nodes.in.play)
    }
    
    tt=tvec[length(tvec)]
  }
  
  ##cleanup
  ixs <- which(!colSums(is.na(HGTs))==nrow(HGTs))
  HGTs <- HGTs[,ixs]
  # if (length(ixs)>0){
  #   rownames(HGTs) <- c('Receiver node','Donor node','time')
  # }
  
  ixx <- which(!colSums(is.na(X))==nrow(X))
  X <- X[,ixx]
  
  return(list('X'=X,'HGTs'=HGTs,'time'=tvec,'tree'=tree))
}