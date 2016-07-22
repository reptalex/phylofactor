#' Colors tree edges based on taxanomic ID
#'
#' @export
#' @param tree phylogenetic tree
#' @param taxonomy list whose first column is OTUs and second column is green-genes style taxonomy, e.g. k__Bacteria; p__Bacteroidetes
#' @param level String indicating taxonomic level for coloring. Must be in set {k,p,c,o,f,g,s}
#' @param outputlegend Logical indicating whether or not to output the list of taxa and their corresponding colors for manual generation of a legend
#' @param colorfcn Color Palette for taxa. Default is rainbow.
#' @param legend Logical whether or not to include legend. Most trees produce pretty ugly legends due to large numbers of taxa, so default is \code{FALSE}
#' @param scramble Logical - whether or not to scramble the order of taxa-color associations by a random draw. Useful for re-drawing colors in case closely related taxa are difficult to distinguish because of adjacent colors in colorfcn.
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
ColorTaxa <- function(tree,taxonomy,level='p',outputlegend=F,colorfcn=NULL,legend=FALSE,scramble=F,...){
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
  
  # Now we clean up the Edgelist a bit to remove redundant edges. 
  atms <- bins(G=Edgelist,set=1:Nedge(tree))
  keepers <- which(Edgelist %in% atms)
  fixers <- setdiff(1:nT,keepers)
  if (length(keepers)==0){stop('Consult Alex Washburne on a funky plot - you encountered No-Mans land')}
  # For each element of "fixers", we need to 
  for (nn in 1:length(fixers)){
    Edgelist[[fixers[nn]]] <- setdiff(Edgelist[[fixers[nn]]],unlist(Edgelist[setdiff(1:nT,fixers[nn])]))
  }

  #Then we assign colors to each taxon
  if (is.null(colorfcn)){
    colorfcn <- rainbow
  }
  colors <- colorfcn(nT)
  if (scramble){
    colors <- colors[sample(nT)]
  }
  
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
