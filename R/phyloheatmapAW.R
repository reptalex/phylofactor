#Ghetto, AW version of phylo.heatmap - optionally does not include tip labels, but does not appropriately re-scale margins.
phylo.heatmapAW <- function(tree,Y,tipLabels=TRUE,...){
  if (all(colnames(Y) %in% tree$tip.label) && all(tree$tip.label %in% colnames(Y)) == FALSE){
    if (all(rownames(Y) %in% tree$tip.label) && all (tree$tip.label %in% rownames(Y)) == TRUE){
      Y <- t(Y)
      } else {
    stop('cannot match rows or columns of Y to tiplabels of tree')}
  edgs <- tree$edge
  n <- length(tree$tip.label)

  #extract order in which taxa show up on tree, from bottom to top
  ord <- match(edgs[,2],1:n)
  ord <- ord[is.na(ord)==F]


  if (n==dim(Y)[1]){Y <- t(Y)}
  m <- dim(Y)[1]

  layout(matrix(c(1,1,2),1,3))
  ### image of data ###
  image(1:m,1:n,Y[,ord],main='Data',xlab="sample",ylab='OTU')

  ### left-facing phylogeny with colored tip-labels corresponding to stacked plot above ###
  par(mar=c(5, 4, 4, 2)+.1)
  pp <- plot.phylo(tree,type='p', use.edge.length=F,direction='l',main='Community Phylogeny',show.tip.label=FALSE,...)
  if (tipLabels==TRUE){
    tiplabels(text= tree$tip.label,bg='white',col = 'black',cex=3)
  }
  ###########################
}
