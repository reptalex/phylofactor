#' Colors tree edges based on taxanomic ID
#'
#' @export
#' @param tree phylogenetic tree
#' @param taxonomy list whose first column is OTUs and second column is green-genes style taxonomy, e.g. k__Bacteria; p__Bacteroidetes
#' @param level String indicating taxonomic level for coloring. Must be in set {k,p,c,o,f,g,s}
#' @param ... additional arguments for plot.phylo
#' @examples
#' data(FTmicrobiome)
#' tree <- FTmicrobiome$PF$tree
#' taxonomy <- FTmicrobiome$taxonomy
#'
#' ColorTaxa(tree,taxonomy,level='p',legend=T,show.tip.label=FALSE,type='unrooted')
#'
#'
#'
#' L <- ColorTaxa(tree,taxonomy,level='c',legend=F,show.tip.label=FALSE,type='unrooted',outputlegend=T)
#'
#' lims <- par('usr')
#' legend(lims[1],lims[4],legend=L$Taxa,fill=L$colors,cex=.6)
ColorTaxa <- function(tree,taxonomy,level='p',outputlegend=F,legend=FALSE,...){
  if (!level %in% c('k','p','c','o','f','g','s')){stop('unknown level - must be a string in the set {k,p,c,o,f,g,s}')}
  #This function produces a tree of the entire community in which the taxa at level "p" are
  #labelled by color.
  if (is.rooted(tree)==F){tree <- root(tree,node=Ntip(tree)+1)}
  #First, we grab the list of taxa
  taxn <- taxonomy[taxonomy[,1] %in% tree$tip.label,]
  Taxa <- listTaxa(taxn,level)

  #Second, for each taxon we grab the OTUids and the unique edges
  nT <- length(Taxa)
  Otus <- lapply(Taxa,grep,x=taxn[,2])
  Edgelist <- vector(mode='list',length=nT)
  for (n in 1:nT){
    Otus[[n]] <- as.character((taxn[Otus[[n]],1]))
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
  ape::plot.phylo(tree,edge.color = edge_colors, edge.width = edge_widths,...)

  if (legend==TRUE){
    lims <- par('usr')
    legend(.99*lims[1],.99*lims[4],gsub(paste(level,'__',sep=''),'',Taxa),fill=colors)
  }
  if (outputlegend){
    output <- list(unlist(Taxa),colors)
    names(output) <- c('Taxa','colors')
    return(output)
  }

}
