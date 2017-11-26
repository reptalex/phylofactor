#' Describes the taxa split at a particular factor summarized by summary.phylofactor
#' @export
#' @param Summary summary.phylofactor object
#' @return list of taxa (with trimmed taxonomy strings) unique to the group & complement for the input summarized factor


TaxaSplit <- function(Summary){
  
  taxG <- as.data.frame(lapply(Summary$group1$IDs,as.character),stringsAsFactors = F)
  taxC <- as.data.frame(lapply(Summary$group2$IDs,as.character),stringsAsFactors = F)

  if (is.null(dim(taxG))){
    nms <- names(taxG)
    taxG <- data.frame('otuIDs'=taxG[1],'TaxaIDs'=taxG[2])
  }
  if (is.null(dim(taxC))){
    nms <- names(taxC)
    taxC <- data.frame('otuIDs'=taxC[1],'TaxaIDs'=taxC[2])
  }
  #### Grab the first unique taxonomic category for the Group and Complement ###

 taxG[,2] <- uniqueTaxa(taxG[,2],taxC[,2])
 taxC[,2] <- uniqueTaxa(taxC[,2],taxG[,2])

 output <- list(taxG,taxC)
 names(output) <- c('Group 1','Group 2')
 return(output)

}
