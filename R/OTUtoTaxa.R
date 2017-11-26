#'Input OTUs and outputs their taxonomic detail, optinally trimmed to their common taxonomic levels.
#'@export
#'@param otus Species IDs found in first column of Taxonomy
#'@param Taxonomy Taxonomy whose first column has OTUs/species IDs and whose second column has greengenes-compatible taxonomic strings
#'@param common.name Logical indicating whether or not to trim the output to the longest common prefix in the taxonomic strings of all otus
#'@param uniques Logical whether or not to trim list to only unique entries
#'@param minimum.tax string indicating the minimum taxonomic detail to include. 
OTUtoTaxa <- function(otus,Taxonomy,common.name=F,uniques=F,minimum.level=Inf){

  if(class(otus)!='character'){otus <- as.character(otus)}
  if (ncol(Taxonomy)>2){
    warning('Taxonomy has more than two columns. Assuming columns 2:ncol(taxonomy) are taxonomic IDs and creating semicolon-delimited taxonomy')
    otus <- Taxonomy[,1]
    tax <- apply(Taxonomy[,2:ncol(taxonomy)],1,paste,sep=';')
    Taxonomy <- data.frame(otus,'taxonomy'=tax)
  }
  if (!class(Taxonomy[,1])=='character'){
    Taxonomy[,1] <- as.character(Taxonomy[,1])
  }
  if (!class(Taxonomy[,2])=='character'){
    Taxonomy[,2] <- as.character(Taxonomy[,2])
  }
  
  notus <- length(otus)
  ind <- match(otus,Taxonomy[,1])
  
  taxa <- listTaxa(Taxonomy[ind,2],minimum.level,common.name,uniques)
  return(taxa)
}
