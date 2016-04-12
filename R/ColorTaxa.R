#' Color tree based on taxonomic level
#' @export
#' @param tree Phylogeny whose tips are contained in the first column of the list, Taxonomy
#' @param Taxonomy list whose first column contains the tip labels of the tree and whose second column contains greengenes taxonomic strings for each OTU.
#' @param level Taxonomic level from set {'k','p','c','o','f','g','s'} default ='p'.
#' @param legend Whether ot not to include a legend of Taxonomic groups at that taxonomic level.
#' @return colored tree labelled by taxa


###################################### ColorTaxa #########################################
ColorTaxa <- function(tree,Taxonomy,level='p',legend=FALSE){
  #This function produces a tree of the entire community in which the taxa at level "p" are
  #labelled by color.
  #First, we grab the list of taxa
  Taxa <- listTaxa(Taxonomy,level)

  #Second, for each taxon we grab the OTUids and the unique edges
  nT <- length(Taxa)
  Otus <- lapply(Taxa,grep,x=Taxonomy[,2])
  Edgelist <- vector(mode='list',length=nT)
  for (n in 1:nT){
    Otus[[n]] <- as.character(Taxonomy[Otus[[n]],1])
    Edgelist[[n]] <- extractEdges(tree,Otus[[n]],type=3)
  }

  #Then we assign colors to each taxon
  colors <- rainbow(nT)
  edge_colors <- rep('black',Nedge(tree))

  for (n in 1:nT){
    edge_colors[Edgelist[[n]]] <- colors[n]
  }
  edge_widths <- rep(2,Nedge(tree))

  # plot the tree
  plot.phylo(tree,use.edge.length = FALSE, show.tip.label = FALSE, type='unrooted',edge.color = edge_colors, edge.width = edge_widths)

  if (legend==TRUE){
    lims <- par('usr')
    legend(.99*lims[1],.99*lims[4],gsub(paste(level,'__',sep=''),'',Taxa),fill=colors)
  }

}
