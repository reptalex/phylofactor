plot.phylofactor <- function(Factor,tree,Data=NULL,bg='white',cex=2,row.names=NULL,clades=1,compare=F,...){
  #returns phylo.heatmap highlighting our factors.
  if (compare==F){
    #makes just one plot

    if (is.null(Data)==T){
      if (is.null(names)){stop('must input rownames for PhyloFactor predicted dataset')}
      Data <- predict.phylofactor(Factor,clades)
      rownames(Data) <- row.names
    }
    phylo.heatmap(tree,clr(Data),...)
    nodelabels(text = as.list(clades),node=Factor$nodes[clades],bg = bg,cex=cex)
  } else {


    #makes two plots for comparison
    par(mfrow=c(2,1))
    if (is.null(Data)==T){
      stop('if compare==T, need to input Data for Comparison')
    }

    PData <-  predict.phylofactor(Factor,clades)
    rownames(PData) <- row.names
    rownames(PData) <- tree$tip.label
    colnames(PData) <- colnames(Data)

    phylo.heatmap(tree,t(clr(t(Data))))
    nodelabels(text = as.list(clades),node=Factor$nodes[clades],bg = bg,cex=cex)
    phylo.heatmap(tree,t(clr(t(PData))))
    nodelabels(text = as.list(clades),node=Factor$nodes[clades],bg = bg,cex=cex)




  }
}
