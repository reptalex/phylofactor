
################################### plot.sig.clades  ############################################
plot.sig.clades <- function(tree,sigOTUs,taxon_site,level='p',legend=FALSE,treetype='unrooted',show.tip.labels=FALSE,use.edge.length=FALSE){
  #this function plots the tree and highlights those edges unique to SigOTUs in the tree,
  # Significant taxa are overlayed on tree labelled by taxonomic level.
  #input tree, sigOTus, taxon_site=taxonomy (mapping from OTU IDs to green genes taxonomy), taxonomic level, whether or not to include a legend, etc.)
  # Uses DJ Bennet's extractEdges function to pull out the edges for sigOTUs. In particular,
  # type = 1: tree if all other taxa besides sigOTUs are dropped
  # type = 2: Edges connecting sigOTUs to terminal node
  # type = 3: edges unique to sigOUTs
  otuids <- tree$tip.label
  if (all(sigOTUs %in% otuids)==FALSE){
    stop('Some SigOTUs are not in the tree')
  }

  Sigedges <- extractEdges(tree,sigOTUs,type=3)

  # plot the tree
  Taxa <- listTaxa(taxon_site,level)

  #Second, for each taxon we grab the OTUids and the unique edges
  nT <- length(Taxa)
  Otus <- lapply(Taxa,grep,x=taxon_site[,2])
  Edgelist <- vector(mode='list',length=nT)
  for (n in 1:nT){
    Otus[[n]] <- as.character(taxon_site[Otus[[n]],1])
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
