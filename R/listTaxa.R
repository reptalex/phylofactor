#' Lists all unique taxonomic names in a Taxonomy at a given level
#' @export
#' @param taxonomy Taxonomy. First row is species IDs, second row is their taxonomy. This function assumes greengenes-compatible taxonomic strings
#' @param level Taxonomic level from set {'k','p','c','o','f','g'}
#' @return list of taxonomic strings for all taxa up to speficied level.
#' @examples
#' data(FTmicrobiome)
#' listTaxa(FTmicrobiome$taxonomy,level='p')

listTaxa <- function(taxonomy,level='p'){
  #This function inputs a taxonomy table whose first column is OTU ids and whose second column is
  #greengenes taxonomy. Output will be a list of the uniqiue taxa at taxonomic level= "level" in the table.
  # level can take arguments in {k,p,c,o,f,g,s}


  pattern1 <- paste(level,'_',sep='')

  if (is.null(dim(taxonomy))==F){
    ntaxa <- dim(taxonomy)[1]
    Taxa <- vector(mode='list',length=ntaxa)

    for (tt in 1:ntaxa){

      s <- toString(taxonomy[tt,2])

      strt <- gregexpr(pattern=pattern1,s)[[1]][1]
      stp <- gregexpr(pattern=';',s)[[1]]
      stp <- stp[min(which(stp>strt))]

      Taxa[[tt]] <- substr(s,1,stp-1)

    }

  } else {


    ntaxa=1
    s <- toString(taxonomy)
    strt <- gregexpr(pattern=pattern1,s)[[1]][1]
    stp <- gregexpr(pattern=';',s)[[1]]
    stp <- stp[min(which(stp>strt))]

    Taxa <- substr(s,1,stp-1)
}

  output <- unique(Taxa)
  return(output)
}
