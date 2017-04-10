#'Input OTUs and outputs their taxonomic detail, optinally trimmed to their common taxonomic levels.
#'@export
#'@param otus Species IDs found in first column of Taxonomy
#'@param Taxonomy Taxonomy whose first column has OTUs/species IDs and whose second column has greengenes-compatible taxonomic strings
#'@param common.name Logical indicating whether or not to trim the output to the longest common prefix in the taxonomic strings of all otus
#'@param uniques Logical whether or not to trim list to only unique entries
#'@param minimum.tax string indicating the minimum taxonomic detail to include. 
OTUtoTaxa <- function(otus,Taxonomy,common.name=F,uniques=F,minimum.tax='p'){
  if (! minimum.tax %in% c('k','p','c','o','f','g','s')){stop('minimum.tax must be a greengenes taxonomic letter, i.e. k, p, c,...')}

    if(typeof(otus)!='character'){otus <- as.character(otus)}
      notus <- length(otus)
      ind <- match(otus,Taxonomy[,1])
      taxa <- sapply(as.list(Taxonomy[ind,2]),FUN=toString)

    mins <- listTaxa(data.frame(otus,taxa),minimum.tax)
    mins <- unique(mins)
    mins[is.na(mins)] <- 'Unassigned'
  
  if(common.name){
    taxa <- substr(taxa[1], start = 1, stop = Biostrings::lcprefix(taxa[1], taxa[length(taxa)]))
  }
  
  
  if(!any(sapply(mins,FUN=function(a,b) grepl(a,b),b=taxa))){taxa <- NULL}
  
  for (s in unlist(mins)){
    if (!any(grepl(s,taxa))){
      taxa <- c(taxa,s)
    }
  }
  
  if (uniques){
    taxa <- unique(taxa)
  }
  
      return(taxa)
}
