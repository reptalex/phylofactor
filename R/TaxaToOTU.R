TaxaToOTU <- function(Taxa,taxonomy){
  #input taxonomy table and list of Taxa. Output list of OTU IDs for each taxon.
  nT <- length(Taxa)
  Otus <- lapply(Taxa,grep,x=taxonomy[,2])
  for (n in 1:nT){
    Otus[[n]] <- as.character(taxonomy[Otus[[n]],1])
  }
  return(Otus)
}
