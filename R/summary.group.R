summary.group <- function(PF,tree,taxonomy,factor,grp){
  #summarizes the OTUids, taxonomic details, data and predictions for an input group of taxa up to a factor level factor.

  output <- NULL
  otuIDs <- tree$tip.label[grp]
  TaxaIDs <- OTUtoTaxa(otuIDs,taxonomy,common.name=F)
  output$IDs <- data.frame(otuIDs,TaxaIDs)

  output$otuData <- PF$Data[grp,]
  output$PF.prediction <- predict.phylofactor(Factor=PF,factors=1:factor)[grp,]
  colnames(output$PF.prediction) <- colnames(PF$Data)
  rownames(output$PF.prediction) <- rownames(PF$Data[grp,])
  return(output)
}
