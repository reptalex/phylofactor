#' trims tree to tip-labels in dataset.
trimTree <- function(Data,tree){
  if ('phylo' %in% class(tree)){
    if (!all(tree$tip.label %in% rownames(Data))){
      tree <- ape::drop.tip(tree,setdiff(tree$tip.label,rownames(Data)))
    }
  }
  return(tree)
}