#' Lists all unique taxonomic names in a Taxonomy at a given level
#' @export
#' @param s Character vector of taxonomic names. First row is species IDs, second row is their taxonomy. This function assumes greengenes-compatible taxonomic strings
#' @param minimum.level Integer. Minimum number of taxonomic levels to include counting from highest to lowest.
#' @param common.name logical whether or not to return the common name in the list of taxa
#' @param uniques logical whether or not to trim taxonomic names to only unique entries
#' @return list of taxonomic strings for all taxa up to speficied level.
#' @examples
#' data(FTmicrobiome)
#' sapply(FTmicrobiome$taxonomy[,2],as.character) %>%
#'    listTaxa(minimum.level=2,uniques=TRUE)

listTaxa <- function(s,minimum.level=1,common.name=F,uniques=F){
  
  stp <- gregexpr(pattern=';',s) %>% sapply(FUN=function(x) if(minimum.level<=length(x)){x[minimum.level]}else{NA} )
  stp[is.na(stp)] <- nchar(s[is.na(stp)])
  Taxa <- substr(s,1,stp)
  
  if(common.name){
    comsub<-function(x) {
      x<-sort(x)
      d_x<-strsplit(x[c(1,length(x))],"")
      der_com<-match(FALSE,do.call("==",d_x))-1
      ifelse(der_com==0,return(character(0)),return(substr(x[1],1,der_com)))
    }
    Taxa <- comsub(Taxa)
  }

  if (uniques){
    output <- unique(Taxa)
  } else {
    output <- Taxa
  }
  return(output)
}
