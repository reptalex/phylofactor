#' Prints phylofactor.summary objects
#' 
#' @export
#' @param output from \code{\link{summary.phylofactor}}
print.phylofactor.summary <- function(s){
  
  info <- s$info
  tbl.str <- paste(capture.output(print.data.frame(s$group.summary)),collapse='\n')
  method <- info$method
  if (!is.null(info$algorithm)){
    algorithm <- paste('
Algorithm                 : ',info$algorithm,sep='')
  } else{
    algorithm <- NULL
  }
  if (!is.null(info$choice)){
    choice <- paste('
Choice                    : ',info$choice,sep='')
  } else {
    choice <- NULL
  }
  
  if (!is.null(info$formula)){
    formula <- paste('
Formula                   : ',info$formula,sep='')
  } else {
    formula <- NULL
  }
  
  nedges <- paste('
Edges Considered          : ',info$nEdges,sep='')
  factor.print <- paste('
Factor                    : ',info$factor,sep='')
  
  factor <- rownames(s$group.summary) %>% strsplit(' ')
  factor <- factor[[1]][2]
  
  n1 <- nrow(s$taxon.tables$Group1)
  tx.tbl1 <- paste(capture.output(print.data.frame(s$taxon.tables$Group1[1:(min(3,n1)),c('Taxon','nSpecies','signal')])),collapse='\n')
  if (n1>3){
    if (n1==4){
      tx.tbl1 <- paste(tx.tbl1,'
                    .................  1 row omitted .................',sep='')
    } else {
      tx.tbl1 <- paste(tx.tbl1,'
                    .................  ',n1-3,' rows omitted .................',sep='')
    }
  }
  n2 <- nrow(s$taxon.tables$Group2)
  tx.tbl2 <- paste(capture.output(print.data.frame(s$taxon.tables$Group2[1:(min(3,nrow(s$taxon.tables$Group2))),c('Taxon','nSpecies','signal')])),collapse='\n')
  if (n2>3){
    if (n2==4){
      tx.tbl2 <- paste(tx.tbl2,'
                    .................  1 row omitted .................',sep='')
    } else {
      tx.tbl2 <- paste(tx.tbl2,'
                    .................  ',n2-3,' rows omitted .................',sep='')
    }
  }
  
  
  ln <- paste('---------------------------------',
              paste(rep('-',nchar(s$info$phylofactor.fcn)),collapse=''),sep='')
  output <- paste('       phylofactor object from function ',s$info$phylofactor.fcn,'
       ',ln,'
Method                    : ',method,algorithm,choice,formula,factor.print,nedges,'

',tbl.str,'
========================================================================
Taxon Tables:
Group1:
',tx.tbl1,'
------------------------------------------------------------------------
Group2:
',tx.tbl2,sep='')
  cat(output)
}
