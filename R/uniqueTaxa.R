#' Find shortest taxonomic prefix
#' @param x list of taxonomic strings whose unique taxonomies we want to know
#' @param y list of taxonomic strings for "subtraction" from x
#' @export
#' @return list of taxonomic strings unique to x, i.e. the shortest taxonomic prefix of x that is unique to x

uniqueTaxa <- function(x,y){
  #returns the first unique taxonomic level for all OTUs in x compared to Y
  lvls <- c('k','p','c','o','f','g')
  nl <- 6
  ntaxa =length(x)
  output <- vector(mode='list',length=ntaxa)
  for (nt in 1:ntaxa){
    tx <- x[nt]

    for (ll in 1:nl){
      ID <- listTaxa(tx,level=lvls[ll])[[1]]
      if (length(grep(ID,y))==0){
        output[[nt]] <- as.factor(ID)
        break
      } else {
        if (ll==nl){
          output[[nt]] <- as.factor(ID)
        }
      }
    }
  }

  return(unlist(output))
}
