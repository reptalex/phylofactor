listTaxa <- function(taxon_site,level='p'){
  #This function inputs a taxonomy table whose first column is OTU ids and whose second column is
  #greengenes taxonomy. Output will be a list of the uniqiue taxa at taxonomic level= "level" in the table.
  # level can take arguments in {k,p,c,o,f,g,s}

  pattern1 <- paste(level,'_',sep='')

  ntaxa <- dim(taxon_site)[1]
  Taxa <- vector(mode='list',length=ntaxa)

  for (tt in 1:ntaxa){

    s <- toString(taxon_site[tt,2])

    strt <- gregexpr(pattern=pattern1,s)[[1]][1]
    stp <- gregexpr(pattern=';',s)[[1]]
    stp <- stp[min(which(stp>strt))]

    Taxa[[tt]] <- substr(s,strt,stp-1)

  }

  output <- unique(Taxa)
  return(output)
}
