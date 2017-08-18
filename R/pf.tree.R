
pf.tree <- function(pf,method='groups',group.map=NA,Grps=NULL,bg.color=NA,bg.alpha=0.1,layout='circular',color.fcn=rainbow,return.ggplot=T,...){
  if (!(is.null(Grps) & is.null(group.map))){
    if (method=='bins'){
      warning('input Grps or group.map will override metohd="bins"')
      method='groups'
    }
  }
  if (method=='groups'){
    if (is.null(Grps)){
      if (is.na(group.map)){
        group.map=data.frame('factor'=1:pf$nfactors,'group'=rep(1,pf$nfactors))
      }
      m <- nrow(group.map)
    } else {
      m <- length(Grps)
    }
  }
  
  n=Ntip(pf$tree)
  
  if (method=='bins'){  
    nd <- numeric(pf$nfactors+1)
    for (i in 1:(pf$nfactors+1)){
      nd[i] <- MRCA(pf$tree,pf$tree$tip.label[pf$bins[[i]]])
    }
    cols <- color.fcn(pf$nfactors+1)
  } else {
    nd <- numeric(m)
    for (i in 1:(m)){
      nd[i] <- MRCA(pf$tree,pf$tree$tip.label[pf$groups[[group.map[i,1]]][[group.map[i,2]]]])
    }
    cols <- color.fcn(m)
  }
  
  ix <- order(nd)                ## this gives us a map of which node to which group/bin
  nd <- sort(nd,decreasing = F)
  
  # gg <- ggtree(pf$tree,...)
  gg <- ggtree(pf$tree,layout=layout)
  if (nd[1]==(n+1) | !is.na(bg.color)){
    gg <- gg+geom_hilight(n+1,fill=bg.color,alpha=bg.alpha)
  }
  
  i=0
  for (ndd in nd){
    i=i+1
    if (!ndd==(n+1)){
      gg <- gg+geom_hilight(ndd,fill=cols[i])
    }
  }
  gg
  
  return(list('map'=ix,'colors'=cols,'ggplot'=gg))
}