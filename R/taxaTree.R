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
