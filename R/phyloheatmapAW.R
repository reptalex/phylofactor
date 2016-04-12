#' Ghetto, AW version of phylo.heatmap
#'@export
#' @param tree phylogeny (ape class)
#' @param Y data matrix. row names must match tip labels of tree
#' @param tipLabels Logical indicating whether or not to include tip labels
#' @param ... optional arguments for plot.phylo(...)
#' @examples
#' set.seed(1)
#' tree <- unroot(rtree(7))

#' X <- as.factor(c(rep(0,5),rep(1,5)))
#' sigClades <- Descendants(tree,c(9,12),type='tips')
#' Data <- matrix(rlnorm(70,meanlog = 8,sdlog = .5),nrow=7)
#' rownames(Data) <- tree$tip.label
#' colnames(Data) <- X
#' Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
#' Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
#' Data <- t(clo(t(Data)))
#'
#' es1 <- extractEdges(tree,tree$tip.label[sigClades[[1]]],type=3)
#' es2 <- extractEdges(tree,tree$tip.label[sigClades[[2]]],type=3)
#' ecs <- rep('black',Nedge(tree))
#' ecs[es1] <- 'red'
#' ecs[es2] <- 'blue'
#' phylo.heatmapAW(tree,Data,tipLabels=F,edge.width=4,edge.color=ecs)
#' nodelabels(c(9,12),c(9,12),cex=2)


#Ghetto, AW version of phylo.heatmap - optionally does not include tip labels, but does not appropriately re-scale margins.
phylo.heatmapAW <- function(tree,Y,tipLabels=TRUE,...){
  if (all(rownames(Y) %in% tree$tip.label) && all (tree$tip.label %in% rownames(Y)) == FALSE){
    if (all(colnames(Y) %in% tree$tip.label) && all(tree$tip.label %in% colnames(Y)) == TRUE){
      Y <- t(Y)
      } else {
    stop('cannot match rows or columns of Y to tiplabels of tree')}}

  ord <- match(tree$tip.label,rownames(Y))
  m <- dim(Y)[1]

  layout(matrix(c(1,1,2),1,3))
  ### image of data ###
  image(t(Y[ord,]),main='Data',xlab="sample",ylab='OTU')

  ### left-facing phylogeny with colored tip-labels corresponding to stacked plot above ###
  par(mar=c(5,4,4,1)+1)
  pp <- ape::plot.phylo(tree,type='p', use.edge.length=F,direction='l',main='Community Phylogeny',show.tip.label=FALSE,...)
  if (tipLabels==TRUE){
    ape::tiplabels(text= tree$tip.label,bg='white',col = 'black',cex=3)
  }
  ###########################
}
