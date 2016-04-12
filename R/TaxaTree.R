#' Construct tree for particular taxonomic group
#' @param tree Ape class phylogeny whose tips are contained in the OTUids of the taxonomy
#' @param taxonomy Taxonomy list whose first column is OTUids and whose second column is taxonomic strings
#' @param Taxon the particular Taxonomic string to include in tree.
#' @return tree whose tips are all the elements in the tree with taxonomy containing the string \code{Taxon}
###################################### TaxaTree ###########################################
TaxaTree <- function(tree,taxonomy,Taxon){
  #This function returns the sub-tree formed by including only those taxa in "Taxon", a greengenes identifier,
  # e.g. p__Firmicutes.
  # Must include taxonomy which maps otuids in tree to the taxonomic labels in Taxon
  OTUs <- as.character(taxonomy[grep(Taxon,taxonomy[,2]),1])
  return(drop.tip(tree,tip = setdiff(tree$tip.label,OTUs)))
}
