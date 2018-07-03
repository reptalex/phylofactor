#' Find shortest unique prefix in two taxonomic lists x and y
#' @param x list of taxonomic strings whose unique taxonomies we want to know
#' @param y list of taxonomic strings for "subtraction" from x
#' @export
#' @return list of taxonomic strings unique to x, i.e. the shortest taxonomic prefix of x that is unique to x

uniqueTaxa <- function(x,y){
  #returns the first unique taxonomic level for all OTUs in x compared to Y

  ntaxa =length(x)
  output <- vector(mode='list',length=ntaxa)
  for (nt in 1:ntaxa){
    tx <- x[nt]
    gg <- gregexpr(';',tx)[[1]]
    nl <- length(gg)
    for (ll in 1:nl){
      # ID <- listTaxa(tx,level=lvls[ll])[[1]]
      ID <- substr(tx,1,gg[ll])
      if (length(grep(ID,y))==0){
        output[[nt]] <- ID
        break
      } else {
        if (ll==nl){
          output[[nt]] <- tx
        }
      }
    }
  }
  
  ## remove final semicolons
  output <- unlist(output)
  ix <- base::endsWith(output,';') & !is.na(output)
  output[ix] <- substr(output[ix],1,nchar(output[ix])-1)
  
  return(output)
}
