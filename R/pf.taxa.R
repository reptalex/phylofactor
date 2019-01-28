#' summarizes taxonomic composition at given factor
#' 
#' @export
#' @param pf phylofactor object
#' @param taxonomy taxonomy whose first column is tree tip labels in phylofactor object and whose second column is semicolon-delimited taxonomy
#' @param factor positive integer. Which factor to summarize.
#' @examples
#' library(phylofactor)
#' data(FTmicrobiome)
#' pf <- FTmicrobiome$PF
#' tx <- FTmicrobiome$taxonomy
#' 
#' pf.taxa(pf,tx,1)
pf.taxa <- function(pf,taxonomy,factor=1){
  if ('data.table' %in% class(taxonomy)){
    taxonomy <- as.data.frame(taxonomy)
  }
  g1 <- pf$tree$tip.label[pf$groups[[factor]][[1]]]
  g2 <- pf$tree$tip.label[pf$groups[[factor]][[2]]]
                            
  t1 <- as.character(taxonomy[match(g1,taxonomy[,1]),2])
  t2 <- as.character(taxonomy[match(g2,taxonomy[,1]),2])
  
  if (all(is.na(t1)) | all(is.na(t2))){
    stop('could not match tree tip-labels of groups to first column of taxonomy')
  }
  
  output <- NULL
  output$group1 <- uniqueTaxa(t1,t2) %>% unique
  output$group2 <- uniqueTaxa(t2,t1) %>% unique
  return(output)
}
