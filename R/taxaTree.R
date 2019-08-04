#' make tree from vector of standardized taxonomies
#' 
#' @param taxa character vector of semicolon-delimited taxonomic names
#' @export
#' @examples 
#' taxa <- c('Prokarya; Bacteria; Proteobacteria; rando1',
#'           'Prokarya; Bacteria; Proteobacteria; rando2',
#'           'Eukarya; Animalia; Vertebrata; cat',
#'           'Eukarya; Animalia; Vertebrata; dog',
#'           'Eukarya; Fungi; Ascomycota; yeast')
#' tree <- taxaTree(taxa)
#' plot(tree)
#' 
#' 
#' library(phylofactor)
#' library(ggplot2)
#' library(ggpubr)
#' data("FTmicrobiome")
#' OTUTable <- FTmicrobiome$OTUTable
#' taxonomy <- FTmicrobiome$taxonomy
#' MetaData <- FTmicrobiome$X
#' 
#' ### need unique taxonomic identifiers for duplicate species
#' taxa <- taxonomy$taxonomy
#' duplicate.taxa <- which(duplicated(taxa))
#' taxa[duplicate.taxa] <- paste(taxa[duplicate.taxa],
#'                               duplicate.taxa,sep='_')
#' 
#' ### With that, we can make a taxonomy tree!
#' tree <- taxaTree(taxa)  ### this tree contains polytomies
#' 
#' ### We need to rename rows of our OTUTable by our unique taxonomy
#' rownames(OTUTable) <- taxa[match(rownames(OTUTable),taxonomy$OTU_ID)]
#' taxonomy$OTU_ID <- taxa ### our unique taxonomic names can be our OTU_IDs
#' 
#' ### phylofactorization
#' pf <- PhyloFactor(OTUTable,tree,MetaData,nfactors=2,ncores=7)
#' 
#' ### Summary and Plotting 
#' 
#' tree.plot <- pf.tree(pf,top.layer = F)
#' clade.colors <- tree.plot$legend$colors
#' taxa1 <- pf.taxa(pf,taxonomy,1)[[1]] ## longest unique taxonomic prefix for our first factor
#' taxa2 <- pf.taxa(pf,taxonomy,2)[[1]]
#' Summary.Table <- data.table('Body_Site'=rep(pf$X,2),
#'                             'ILR_Abundance'=c(pf$models[[1]]$y,
#'                                               pf$models[[2]]$y),
#'                             'taxa'=rep(c(taxa1,taxa2),
#'                                        each=ncol(OTUTable)))
#' abundance.plot <- ggplot(Summary.Table,aes(Body_Site,ILR_Abundance,color=taxa))+
#'   geom_boxplot()+
#'   geom_jitter(cex=2)+
#'   facet_wrap(.~taxa,nrow=2)+
#'   scale_x_discrete('Body Site')+
#'   scale_y_continuous('ILR Abundance')+
#'   scale_color_manual(values = clade.colors)+
#'   theme(legend.position='none')
#' 
#' ggarrange(tree.plot$ggplot,abundance.plot,ncol=2,widths = c(1,2))
taxaTree <- function(taxa){
  taxa.list <- sapply(taxa,strsplit,';')
  
  if (length(unique(sapply(taxa.list,length)))>1){
    stop('taxa vary in length - please standardize labels for comparison')
  }
  
  taxa.similarity <- function(t1,t2){
    d=0
    for (i in 1:min(length(t1),length(t2))){
      ix <- as.numeric(t1[[i]]==t2[[i]])
      d=d+ix
      if (ix==0){
        break
      }
    }
    return(d)
  }
  
  m <- max(sapply(taxa.list,length))
  D <- matrix(m,
              nrow=length(taxa.list),
              ncol=length(taxa.list))
  for (i in 1:(length(taxa.list)-1)){
    for (j in (i+1):length(taxa.list)){
      n <- taxa.similarity(taxa.list[[i]],taxa.list[[j]])
      D[i,j] <- n
      D[j,i] <- n
    }
  }
  
  D <- m-D
  tree <- ape::as.phylo(stats::hclust(stats::as.dist(D)))
  tree$tip.label <- taxa
  tree <- ape::di2multi(tree)
  return(tree)
}
