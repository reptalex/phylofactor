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
