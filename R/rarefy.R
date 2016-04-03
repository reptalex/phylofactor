rarefy <- function(OTUTable,to='min',method='multinom'){
  
  if (to=='min'){
    lvl <- min(colSums(OTUTable))
  } else if (to=='max'){
    lvl <- max(colSums(OTUTable))
  } else {
    if (is.numeric(to)==F){stop('unknown input, "to"')}
  }
  
  if (method=='multinom'){
    output <- apply(OTUTable,MARGIN=2,FUN=function(x,lvl){return(rmultinom(1,lvl,x))},lvl=lvl)
  } else if (method=='nbinom') {
    
  }
  
  colnames(output) <- colnames(OTUTable)
  rownames(output) <- rownames(OTUTable)
  return(output)
}