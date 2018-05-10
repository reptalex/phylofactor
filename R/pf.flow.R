#' Generates flow of phylofactor object
#' 
#' @param PF phylofactor class object
#' @param factor integer to which factor the flow should detail

pf.flow <- function(PF,factor=PF$nfactors){
  b <- PF$Var$basis[,1:factor,drop=F]
  bs <- bins(b)
  output <- NULL
  nt <- length(bs)
  dmat <- matrix(factor+1,ncol=factor+1,nrow=factor+1)
  diag(dmat) <- 0
  for (ff in rev(1:factor)){
    v <- b[,ff]
    ix <- which(v!=0)
    ins <- which(sapply(bs,FUN=function(a) all(a %in% ix)))
    ins <- combn(ins,2,simplify=FALSE)
    for (i in 1:length(ins)){
      dmat[ins[[i]][1],ins[[i]][2]]=dmat[ins[[i]][1],ins[[i]][2]]-1
      dmat[ins[[i]][2],ins[[i]][1]]=dmat[ins[[i]][2],ins[[i]][1]]-1
    }
  }
  colnames(dmat) <- output$tip.label
  rownames(dmat) <- output$tip.label
  output <- ape::as.phylo(stats::hclust(dist(dmat)))
  output$tip.label <- sapply(as.list(output$tip.label),FUN=function(a) paste('Bin',a))
  
  ix <- sapply(as.list(1:nt),FUN=function(a,b) which(b==a),b=output$edge[,2]) %>% order
  
  output$order <- output$tip.label[ix]
  
  return(output)
}