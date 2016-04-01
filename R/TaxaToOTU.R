TaxaToOTU <- function(taxon_table,Taxa){
  #input taxonomy table and list of Taxa. Output list of OTU IDs for each taxon.
  nT <- length(Taxa)
  Otus <- lapply(Taxa,grep,x=taxon_site[,2])
  for (n in 1:nT){
    Otus[[n]] <- as.character(taxon_site[Otus[[n]],1])
  }
  return(Otus)
}
