plot.phylofactor <- function(Factor,tree,Data=NULL,nclades=1,...){
  #returns phylo.heatmap highlighting our factors. 
  
  if (is.null(Data)==T){
    Data <- predict.phylofactor(Factor)
  }
  phylo.heatmap(tree,clr(Data),...)
  nodelabels(text = as.list(1:nclades),node=Factor$nodes[1:nclades])
}