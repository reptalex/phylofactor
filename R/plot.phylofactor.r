#' Plots phylogeny and the data, highlighting nodes pulled out by PhyloPF object and optionally comparing to the PhyloPF prediction.
#'
#'  @param PF PhyloFactor object
#'  @param tree Phylogeny used in generating PhyloFactor object. Default will scan PF for a tree.
#'  @param Data Data used to generate PhyloFactor object. If Null, default will plot the phylofactor prediction of the dataset.
#'  @param bg Background color for node labels, default 'white'
#'  @param cex Size of node labels, default 2.
#'  @param clades Clades, also referred to as the splits or factors, of phylofactor object to be labelled on tree and used in prediction.
#'  @param compare Whether or not to compare the Data with phylofactor prediction. If compare = TRUE, two phylo-heatmaps will be produced for direct comparison.
#'  @param ... additional arguments for phylo.heatmap
#'  @example
#'  set.seed(1)
#'  tree <- rtree(7)
#'  X <- as.factor(c(rep(0,5),rep(1,5)))
#'  sigClade <- Descendants(tree,10,type='tips)
#'  Data <- matrix(rlnorm(70),nrow=7)
#'  rownames(Data) <- tree$tip.label
#'  colnames(Data) <- X
#'  Data[sigClade,X==1] <- Data[sigClade,X==1]*8
#'  Data <- t(clo(t(Data)))
#'
#'  PF <- PhyloFactor(Data,tree,nclades=1)
#'
#'  plot.phylofactor(PF,compare=T)
#'
#'
#'  Data[sigClade,]

plot.phyloPF <- function(PF,tree=NULL,Data=NULL,bg='white',cex=2,clades=1,compare=F,...){
  #returns phylo.heatmap highlighting our PFs.
  if(is.null(tree)){
    if(is.null(PF$tree)){stop('Input Phylofactor object does not contain tree - must input tree')}
    tree <- PF$tree
  }

  if (compare==F){
    #makes just one plot

    if (is.null(Data)==T){
      if (is.null(names)){stop('must input rownames for PhyloPF predicted dataset')}
      row.names=tree$tip.label
      PData <- predict.phyloPF(PF,clades)
      rownames(PData) <- row.names
      colnames(PData) <- colnames(PData)
      phylo.heatmap(tree,t(clr(t(PData))),...)
    } else {
      phylo.heatmap(tree,t(clr(t(Data))),...)
    }
    nodelabels(text = as.list(clades),node=PF$nodes[clades],bg = bg,cex=cex)
  } else {


    #makes two plots for comparison
    par(mfrow=c(2,1))
    if (is.null(Data)==T){
      stop('if compare==T, need to input Data for Comparison')
    }
    row.names=rownames(Data)
    PData <-  predict.phyloPF(PF,clades)
    rownames(PData) <- row.names
    rownames(PData) <- tree$tip.label
    colnames(PData) <- colnames(Data)

    phylo.heatmap(tree,t(clr(t(Data))))
    nodelabels(text = as.list(clades),node=PF$nodes[clades],bg = bg,cex=cex)
    phylo.heatmap(tree,t(clr(t(PData))))
    nodelabels(text = as.list(clades),node=PF$nodes[clades],bg = bg,cex=cex)




  }
}
