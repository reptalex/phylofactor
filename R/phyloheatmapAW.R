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
  pp <- plot.phylo(tree,type='p', use.edge.length=F,direction='l',main='Community Phylogeny',show.tip.label=FALSE,...)
  if (tipLabels==TRUE){
    tiplabels(text= tree$tip.label,bg='white',col = 'black',cex=3)
  }
  ###########################
}
