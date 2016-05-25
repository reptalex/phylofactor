

atomTaxa <- function(PF,nfactors,taxonomy,common.name=F,uniques=T){
  atms <- atoms(PF$basis[,1:nfactors,drop=F])
  otus <- lapply(atms,FUN = function(x,names) names[x],names=PF$tree$tip.label)
  taxa <- lapply(otus,FUN=OTUtoTaxa,Taxonomy=taxonomy,common.name=common.name,uniques=uniques)
  taxa2 <- taxa
  
  if (uniques){
    for (nn in 1:(nfactors+1)){
      taxa2[[nn]] <- unique(uniqueTaxa(taxa[[nn]],unlist(taxa[setdiff(1:length(taxa),nn)])))
    }
    taxa=taxa2
  }
  
  return(taxa)
}