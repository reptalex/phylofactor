#' Visualize results form simHGT
#' 
#' @export
#' @param sim Simulation output from \code{\link{simHGT}}
#' @param main title of plot
#' @param cex cex applied across plot, points and lines.
#' @param lwd line width
#' @param pch point choice, see \code{\link{plot}}
#' @examples
#' library(phylofactor)
#' 
#' set.seed(5)
#' tree <- rtree(5)
#' sim <- simHGT(tree,0.3)
#' par(mfrow=c(2,1))
#' plot.phylo(tree)
#' HGTplot(sim,lwd=2,main='HGT Dynamics')

HGTplot <- function(sim,main='Trait Evolution',cex=1.5,lwd=2,pch=16,...){

X <- sim$X
HGTs <- sim$HGTs
tvec <- sim$time
tree <- sim$tree

hgtline <- function(hgt,tree){
  # plots dashed line with dot indicating HGT event.
  ixr <- hgt[1]  #receiver node
  ixd <- hgt[2]  #donor node
  tt <- hgt[3]   #time
  ixt <- which(tvec==tt)
  
  # get trait values
  if (ixd<=ape::Ntip(tree)){
    lbl <- tree$tip.label[ixd]
  } else {
    lbl <- tree$tip.label[phangorn::Descendants(tree,ixd,'tips')[[1]][1]]
  }
  xd <- X[lbl,ixt[1]] #donors are the same
    
  if (ixr<=ape::Ntip(tree)){
    lbl <- tree$tip.label[ixr]
  } else {
    lbl <- tree$tip.label[phangorn::Descendants(tree,ixr,'tips')[[1]][1]]
  }
  xr <- X[lbl,ixt]
  lines(tvec[ixt],xr,col='red',cex=cex,lwd=lwd)
  points(tvec[ixt[2]],xd,col='black',cex=cex,pch=pch)
}
 

ix <- min(which(apply(X,MARGIN=1,FUN=function(x) sum(!is.na(x)))==ncol(X)))
plot(tvec,X[tree$tip.label[ix],],type='l',ylim=c(min(X,na.rm=T),max(X,na.rm=T)),xlab='time',ylab='Trait Value',main=main,lwd=lwd,...)

for (tip in 1:Ntip(tree)){
  x <- X[tree$tip.label[tip],]
  if (any(is.na(x))){
    ix <- 1:min(which(is.na(x)))
    lines(tvec[ix],x[ix],lwd=lwd)
  } else {
    lines(tvec,x,lwd=lwd)
  }
}

if (length(HGTs)>0){ 
  #We will draw dashed, vertical lines pointing to the HGT events
  if (is.null(dim(HGTs))){#only one HGT
      hgtline(HGTs,tree)
  } else {
    for (n in 1:ncol(HGTs)){
      hgtline(HGTs[,n],tree)
    }
  }
}

}