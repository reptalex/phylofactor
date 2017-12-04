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
#' ColorTaxa(tree,taxonomy,legend=TRUE,show.tip.label=FALSE,type='unrooted')
#'
#' L <- ColorTaxa(tree,taxonomy,level='c',legend=FALSE,show.tip.label=FALSE,
#'                      type='unrooted',outputlegend=TRUE)
#'
#' lims <- par('usr')
#' legend(lims[1],lims[4],legend=L$Taxa,fill=L$colors,cex=.6)
ColorTaxa <- function(tree,taxonomy,minimum.level=2,outputlegend=FALSE,colorfcn=NULL,legend=FALSE,scramble=FALSE,...){

  if (ape::is.rooted(tree)==FALSE){tree <- ape::root(tree,node=ape::Ntip(tree)+1)}
  #First, we grab the list of taxa
  taxn <- taxonomy[taxonomy[,1] %in% tree$tip.label,]
  Taxa <- listTaxa(taxn[,2],minimum.level,uniques = T)

  #Second, for each taxon we grab the OTUids and the unique edges
  nT <- length(Taxa)
  Otus <- lapply(Taxa,grep,x=taxn[,2])
  Edgelist <- vector(mode='list',length=nT)
  for (n in 1:nT){
    Otus[[n]] <- as.character((taxn[Otus[[n]],1]))
    Edgelist[[n]] <- extractEdges(tree,Otus[[n]],type=3)
  }
  
  # Now we clean up the Edgelist a bit to remove redundant edges. 
  atms <- bins(G=Edgelist,set=1:(ape::Nedge(tree)))
  keepers <- which(Edgelist %in% atms)
  fixers <- setdiff(1:nT,keepers)
  if (length(keepers)==0){stop('Consult Alex Washburne on a funky plot - you encountered No-Mans land')}
  # For each element of "fixers", we need to 
  if (length(fixers)>0){
    for (nn in 1:length(fixers)){
      Edgelist[[fixers[nn]]] <- setdiff(Edgelist[[fixers[nn]]],unlist(Edgelist[setdiff(1:nT,fixers[nn])]))
    }
  }

  #Then we assign colors to each taxon
  if (is.null(colorfcn)){
    colorfcn <- rainbow
  }
  colors <- colorfcn(nT)
  if (scramble){
    colors <- colors[sample(nT)]
  }
  
  edge_colors <- rep('black',ape::Nedge(tree))
  for (n in 1:nT){
    edge_colors[Edgelist[[n]]] <- colors[n]
  }
  edge_widths <- rep(2,ape::Nedge(tree))

  # plot the tree
  ape::plot.phylo(tree,edge.color = edge_colors, edge.width = edge_widths,...)

  if (legend==TRUE){
    lims <- par('usr')
    legend(.7*lims[2],.99*lims[4],Taxa,fill=colors)
  }
  if (outputlegend){
    output <- list(unlist(Taxa),colors)
    names(output) <- c('Taxa','colors')
    return(output)
  }

}
