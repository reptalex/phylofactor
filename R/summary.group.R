#' Summarizes OTU IDs, taxonomic detail, data and phylofactor predictions of a given input group of taxa
#' @export
#' @param PF PhyloFactor object
#' @param tree phylogeny (ape class)
#' @param taxonomy First column is OTU ids also found in tree$tip.label and second column is greengenes taxonomic strings
#' @param factor Factor number up to which phylofactor predictions should be used
#' @param grp input vector of taxa to be summarized
#' @return summary object with $IDs, $otuData, $PF.prediction

summary.group <- function(PF,tree,taxonomy,factor,grp){
  #summarizes the OTUids, taxonomic details, data and predictions for an input group of taxa up to a factor level factor.

  output <- NULL
  otuIDs <- tree$tip.label[grp]
  TaxaIDs <- OTUtoTaxa(otuIDs,taxonomy,common.name=F)
  output$IDs <- data.frame(otuIDs,TaxaIDs)

  output$otuData <- PF$Data[grp, ,drop=F]
  output$PF.prediction <- phylofactor.predict(PF,factors=1:factor)[grp, ,drop=F]
  colnames(output$PF.prediction) <- colnames(PF$Data)
  rownames(output$PF.prediction) <- rownames(PF$Data[grp, ,drop=F])
  return(output)
}
