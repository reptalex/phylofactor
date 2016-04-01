###################################### TaxaTree ###########################################
TaxaTree <- function(tree,taxon_site,Taxon){
  #This function returns the sub-tree formed by including only those taxa in "Taxon", a greengenes identifier,
  # e.g. p__Firmicutes.
  # Must include taxon_site which maps otuids in tree to the taxonomic labels in Taxon
  OTUs <- as.character(taxon_site[grep(Taxon,taxon_site[,2]),1])
  return(drop.tip(tree,tip = setdiff(tree$tip.label,OTUs)))
}
