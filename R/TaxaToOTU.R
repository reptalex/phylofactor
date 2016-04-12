#' Quick conversino from taxonomic string to OTU ID
#' @export
#' @param Taxa Taxonomic string, e.g. k__Animalia or 'Animal'
#' @param taxonomy Greengenes format taxonomy whose first column is OTUids and whose second column is taxonomic strings
#' @return list of OTU IDs whose taxonomic string contains "Taxa" input.


TaxaToOTU <- function(Taxa,taxonomy){
  #input taxonomy table and list of Taxa. Output list of OTU IDs for each taxon.
  nT <- length(Taxa)
  Otus <- lapply(Taxa,grep,x=taxonomy[,2])
  for (n in 1:nT){
    Otus[[n]] <- as.character(taxonomy[Otus[[n]],1])
  }
  return(Otus)
}
