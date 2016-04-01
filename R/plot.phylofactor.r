plot.phylofactor <- function(Factor,tree,Data=NULL,clades=1,compare=F,...){
  #returns phylo.heatmap highlighting our factors.
  if (compare==F){
    #makes just one plot

    if (is.null(Data)==T){
      Data <- predict.phylofactor(Factor,clades)
    }
    phylo.heatmap(tree,clr(Data),...)
    nodelabels(text = as.list(clades),node=Factor$nodes[clades])
  } else {


    #makes two plots for comparison
    par(mfrow=c(2,1))
    if (is.null(Data)==T){
      stop('if compare==T, need to input Data for Comparison')
    }

    PData <-  predict.phylofactor(Factor,clades)
    rownames(PData) <- tree$tip.label
    colnames(PData) <- colnames(Data)

    phylo.heatmap(tree,t(clr(t(Data))))
    nodelabels(text = as.list(clades),node=Factor$nodes[clades])
    phylo.heatmap(tree,t(clr(t(PData))))
    nodelabels(text = as.list(clades),node=Factor$nodes[clades])




  }
}
