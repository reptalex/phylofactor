#' Plots tree and highlights otus given by sigOTUs
#'
#' @param tree phylogeny (ape class)
#' @param sigOTUs taxa to highlight, corresponding to tip labels of tree.
#' @param Taxonomy Taxonomy whose first column is OTU ids and whose second column is a greengenes style taxonomic string
#' @param level Taxonomic level from {'k','p','c','o','f','g'} for coloring tree
#' @param legend Logical whether or not to include legend of taxonomic types corresponding to colors in tree
#' @param ... optional arguments for plot.phylo()
#' @return colored plot of phylogeny, highlighting sigOTUs
#' @example
#' data("FTmicrobiome")
#' set.seed(1)
#' tree <- drop.tip(tree,tree$tip.label[sample(Ntip(tree),Ntip(tree)-100)])
#' otus <- tree$tip.label
#' Taxonomy <- FTmicrobiome$taxonomy[which(FTmicrobiome$taxonomy[,1] %in% otus),]
#' txa <- listTaxa(Taxonomy,level='c')
#' txa <- txa[[3]]  #these are the p__Firmicutes; c__Bacilli
#' otus <- TaxaToOTU(txa,Taxonomy)[[1]]
#'
#' plot.sig.clades(tree,otus,Taxonomy,level='p',legend=T)


################################### plot.sig.clades  ############################################
plot.sig.clades <- function(tree,sigOTUs,Taxonomy,level='p',legend=FALSE,...){
  #this function plots the tree and highlights those edges unique to SigOTUs in the tree,
  # Significant taxa are overlayed on tree labelled by taxonomic level.
  #input tree, sigOTus, Taxonomy=taxonomy (mapping from OTU IDs to green genes taxonomy), taxonomic level, whether or not to include a legend, etc.)
  # Uses DJ Bennet's extractEdges function to pull out the edges for sigOTUs.
  otuids <- tree$tip.label
  if (all(sigOTUs %in% otuids)==FALSE){
    stop('Some SigOTUs are not in the tree')
  }

  if (is.list(sigOTUs)==F){
    Sigedges <- extractEdges(tree,sigOTUs,type=3)
  } else {
    Sigedges <- lapply(sigOTUs,FUN=extractEdges,tree=tree,type=3)
  }

  # plot the tree
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
  edge_colors <- rep('red',Nedge(tree))

  for (n in 1:nT){
    edge_colors[Edgelist[[n]]] <- colors[n]
  }
  edge_colors[extractEdges(tree,sigOTUs,type=3)] <- 'black'

  edge_widths <- rep(2,Nedge(tree))
  edge_widths[extractEdges(tree,sigOTUs,type=3)] <- 8

  # plot the tree
  plot.phylo(tree,use.edge.length = FALSE, show.tip.label = FALSE, type='unrooted',edge.color = edge_colors, edge.width = edge_widths)

  if (legend==TRUE){
    lims <- par('usr')
    Taxa <- c(Taxa,'Significant Clades')
    legend(.99*lims[1],.99*lims[4],gsub(paste(level,'__',sep=''),'',Taxa),fill=c(colors,'black'))
  }
}
