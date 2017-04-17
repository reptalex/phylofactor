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
#' set.seed(5)
#' tree <- rtree(5)
#' sim <- simHGT(tree,0.3)
#' par(mfrow=c(2,1))
#' plot.phylo(tree)
#' HGTplot(sim,lwd=2,main='HGT Dynamics')

simHGT <- function(tree,rate,a=0,mu=0,sig2=1){
  # Gillespie simulation of HGT on a brownian motion process
  n <- length(tree$tip.label)                 #Number of extant species
  tipheights <- diag(ape::vcv.phylo(tree))    #distances of species to root
  Tmax <- max(tipheights)                     #maximum distance to root - the maximum time for our sim
  
  nodes <- ape::Ntip(tree)+2:ape::Nnode(tree) #nodes in tree
                                              #groups corresponding to nodes:
  grps <- sapply(as.list(nodes),FUN=function(n,t) phangorn::Descendants(t,n,type='tips'),t=tree)
  names(grps) <- nodes
                                              #distances of nodes to root
  ndheights <- sapply(as.list(nodes),FUN=function(n,t) phytools::nodeheight(t,n),t=tree)
  names(ndheights) <- nodes
  
  
  approx.hgt <- ceiling(rate*sum(tree$edge))  #approximate number of HGT events we will see
  HGTs <- matrix(NA,nrow=3,ncol=2*approx.hgt) #Matrix documenting HGT events: receiver/donor nodes and timing
  
                                              #State matrix
  X <- matrix(NA,ncol=(2*n-1)+2*approx.hgt,nrow=n)
  rownames(X) <- tree$tip.label
  X[,1] <- a                                  #Initialize at ancestral state
  binder <- matrix(NA,ncol=1,nrow=nrow(X))    #Used to grow X as needed
  
  
  ############################## initialization ###########################
  tt <- 0     #stopwatch
  event <- 1  #event counter
  hgt <- 0    #hgt event counter
  tvec <- 0   #vector of event times, will be equal to ncol(X).
              #nodes.in.play: the extant nodes at time tt, 
              #before which diffusion can occur and between which HGT may occur 
  nodes.in.play <- phangorn::Descendants(tree,n+1,type='children') 
  nn <- length(nodes.in.play)  # Number of edges currently evolving
  nt <- length(tipheights)     # Number of tips still in play - extinct tips no must not diffuse.
  
  while (tt<Tmax){
    
    
    ## Find distance to next tip or node
    possibilities <- c(tipheights,ndheights)
    ix.min <- which.min(possibilities)
    dtmax <- possibilities[ix.min]-tt  #maximum dt from now until we must update nodes.in.play and tipheights
    
    ## Number of HGTs that happen between now and next tip or node
    n.hgts <- rpois(1,lambda=rate*dtmax*nn)
    if (n.hgts>0 & length(nodes.in.play)==1){n.hgts=0}
    
    if (n.hgts>0){
      ######### HGT events happened ##########
      
      ##### timing of events
      ts <- runif(n.hgts,min=tt,max=tt+dtmax) %>% sort
      ### Note: since there is a constant number of extant lineages over our window, the timing of 
      ### poisson HGT events will be uniformly distributed over [tt,tt+dtmax]
      
      for (i in 1:n.hgts){
        event=event+1   ##Event: drift before HGT
        hgt=hgt+1       ##update hgt for HGTs
        
        
        ############### Diffusion between HGT events ####################
        ##### Time of diffusion: DT
        if (i==1){
          DT <- ts[i]-tt
        } else {
          DT <- ts[i]-ts[i-1]
        }
        tvec <- c(tvec,ts[i])
        
        ####### Update X ######################################
        
        if (event>ncol(X)){
          X <- cbind(X,binder)
        }
        
        ####### Diffusion over DT
        dx <- rnorm(nn,mean=mu*DT,sd = sqrt(sig2*DT)) 
        
        for (ll in 1:length(nodes.in.play)){
          nde <- nodes.in.play[ll]
          if (nde>n){ #node is an internal node
            ##### get indexes of species being modified
            ixs <- tree$tip.label[grps[toString(nde)][[1]]]
          } else {    #node is a tip
            ixs <- tree$tip.label[nde]
          }
          X[ixs,event] <- X[ixs,event-1]+dx[ll]
        }

        ########################################################
        
        event=event+1 ### Event: HGT
        tvec <- c(tvec,ts[i]) ##note: here, tvec[event]=tvec[event-1]; this allows graphical illustration of HGT
        ########### HGT ########################################
        ### Randomly sample receiver and donor node
        ### and update HGTs to include receiver/donor node and time
        if (hgt>ncol(HGTs)){
          hb <- c(sample(nodes.in.play,2,replace=F),ts[i])
          HGTs <- cbind(HGTs,matrix(hb,ncol=1))
        } else {
          HGTs[1:2,hgt] <- sample(nodes.in.play,2,replace=F)
          HGTs[3,hgt] <- ts[i]
        }
        
        ####### Get indexes of receiver and donors ############
        if (HGTs[1,hgt]<=n){ # receiving node is a tip
          ixr <- tree$tip.label[HGTs[1,hgt]]
        } else {             # receiving node is internal
          ixr <- tree$tip.label[grps[toString(HGTs[1,hgt])][[1]]] #receiving node
        }
        if (HGTs[2,hgt]<=n){ # donor node is a tip
          ixd <- tree$tip.label[HGTs[2,hgt]]
        } else {             # donor node is internal
          ixd <- tree$tip.label[grps[toString(HGTs[2,hgt])][[1]]] #receiving node
        } 
        
        ################## Update X and tvec ##################
        if (event>ncol(X)){
          X <- cbind(X,binder)
        }
        
        ### Receiver      Donor 
        X[ixr,event] <- X[ixd[1],event-1]  ### replace receiver traits with donor traits
        ### All  non-receiving states remain unchanged
        ixn <- setdiff(rownames(X),ixr)
        X[ixn,event] <- X[ixn,event-1]
      }  ######### END HGT LOOP ########
      
      
      ## simulate the last bit of drift until our next split or terminal node
      event=event+1   ### Event: drift after HGT
      if (event>ncol(X)){
        X <- cbind(X,binder)
      }
      tvec <- c(tvec,tt+dtmax)
      DT <- (tt+dtmax)-ts[i] #final DT
      dx <- rnorm(nn,mean=mu*DT,sd = sqrt(sig2*DT))
      for (ll in 1:length(nodes.in.play)){
        nde <- nodes.in.play[ll]
        if (nde>n){
          ixs <- tree$tip.label[grps[toString(nde)][[1]]]
        } else {
          ixs <- tree$tip.label[nde]
        }
        X[ixs,event] <- X[ixs,event-1]+dx[ll]
      }
      
    } else { ### no HGT
      
      event=event+1 ### Event: Diffusion over dtmax
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
    
    
    ######### INCOMPLETE - Extinct lineages still change. Need to fix that.
    if (ix.min <= nt){ #we finish at a tip. Need to remove it from the pool
      tipheights <- tipheights[setdiff(1:nt,ix.min)]
        # setdiff(tipheights,tipheights[ix.min])
      ## remove tip from nodes.in.play
      ix <- match(names(ix.min),tree$tip.label)
      nodes.in.play <- setdiff(nodes.in.play,ix)
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
  if (length(ixs)>1){
    rownames(HGTs) <- c('Receiver node','Donor node','time')
  } else if (length(ixs)==1){
    names(HGTs) <- c('Receiver node','Donor node','time')
  } 

  ixx <- which(!colSums(is.na(X))==nrow(X))
  X <- X[,ixx]
  
  states <- apply(X,MARGIN=1,FUN=function(x) x[max(which(!is.na(x)))])
  
  return(list('X'=X,'states'=states,'HGTs'=HGTs,'time'=tvec,'tree'=tree))
}