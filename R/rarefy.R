#' Rarefies and OTU table to a given level
#' @export
#' @param OTUTable Table of species abundances/counts whose rows are species and whose columns are samples
#' @param to Number of counts post-rarefaction, from 'min','max', or an integer. Default is 'min', minimum number of counts of all samples in the table
#' @param method Method for randomly drawing counts. Currently only supports "multinom", but coming soon: 'nbinom' and 'PoisLN'
#' @return Table of abundances with colSums(Table) = to number specified.

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
