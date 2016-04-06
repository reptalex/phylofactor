TaxaSplit <- function(NodeSummary){

  #
  taxG <- NodeSummary$group$IDs
  taxC <- NodeSummary$complement$IDs

  #### Grab the first unique taxonomic category for the Group and Complement ###

 taxG[,2] <- uniqueTaxa(taxG[,2],taxC[,2])
 taxC[,2] <- uniqueTaxa(taxC[,2],taxG[,2])

 output <- list(taxG,taxC)
 names(output) <- c('Group','Complement')
 return(output)

}
