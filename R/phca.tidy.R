#' clean description of phca object
#' 
#' @export
#' @param phca \code{PhyCA} object
#' @param Taxonomy Greengenes taxonomy, first column are otus and second column are taxonomic strings of format "p__Phylum; c__Class; o__Order..."
#' @param ncomponents Integer number of PhyCA components to tidy, from \code{1:ncomponents}
#' @param common.name Logical, whether or not to trim output taxonomic summary to the longest common prefix. Default \code{FALSE}
#' @param uniques Logical, whether or not to trim taxonomic output to unique taxonomic names. Default \code{TRUE}
#' @param getEdges Logical, whether or not to get edges corresponding to the splits. Default \code{FALSE}
#' @param ncores Integer number of cores to use if getting edges - does not affect runtime if \code{getEdges=FALSE}
#' @examples
#' 
#' library(phylofactor)
#' data("FTmicrobiome")
#' 
#' Data <- FTmicrobiome$PF$Data
#' tree <- FTmicrobiome$PF$tree
#' taxonomy <- FTmicrobiome$taxonomy
#' 
#' phca <- PhyCA(Data,tree,ncomponents = 3,ncores=3,output.edges=F)
#' 
#' #the standard
#' td <- phyca.tidy(phca,Taxonomy)
#' 
#' #or, for a simpler view
#' td <- phyca.tidy(phca,taxonomy,taxa.split=T)
#' 
#' sum(FTmicrobiome$PF$basis[,1:3]-phca$basis[,1:3])  #the first three phylofactors here correspond to the first three ILR-phylogenetic Components. 

phyca.tidy <- function(phca,Taxonomy,ncomponents=NULL,taxa.split=F,common.name=F,uniques=T,getEdges=F,ncores=NULL,...){
  
  if (is.null(ncomponents)){
    ncomponents <- ncol(phca$basis)
  }
  OTUs <- rownames(phca$Data)
  output <- NULL
  output$taxa <- vector(mode='list',length=ncomponents)
  
  output$summary <- data.frame('Group 1' = numeric(ncomponents),'Group 2'=numeric(ncomponents),'Percent Variance'=phca$PercentVariance[1:ncomponents])
  rownames(output$summary) <- lapply(as.list(1:ncomponents),FUN=function(s) paste('phyComp ',toString(s),':',sep=''))
  if (getEdges){
    output$edges <- vector(mode='list',length=ncomponents)
  }

  
  for (nn in 1:ncomponents){
    x <- phca$basis[,nn]
    grp <- list(which(x>0),which(x<0))
    grp <- getLabelledGrp(tree=phca$tree,Groups=grp)
    otus <- lapply(grp,FUN = function(g,otus) otus[g],otus=OTUs) 
    
    if (!taxa.split){
      taxa <- lapply(otus,FUN=function(g,Taxonomy,common.name,uniques,...) OTUtoTaxa(g,Taxonomy=Taxonomy,common.name,uniques,...),common.name=common.name,uniques=uniques,Taxonomy=Taxonomy)
      output$taxa[[nn]] <- taxa
    } else {
      taxa <- lapply(otus,FUN=function(g,Taxonomy,common.name,uniques,...) OTUtoTaxa(g,Taxonomy=Taxonomy,common.name,uniques,...),common.name=F,uniques=T,Taxonomy=Taxonomy)
      output$taxa[[nn]] <- vector(mode='list',length=2)
      output$taxa[[nn]][[1]] <- unique(uniqueTaxa(taxa[[1]],taxa[[2]]))
      output$taxa[[nn]][[2]] <- unique(uniqueTaxa(taxa[[2]],taxa[[1]]))
    }
    output$summary[nn,1:2] <- sapply(as.list(names(grp)),FUN=function(a) paste('--(',a,')--',sep=''))
    names(output$taxa[[nn]]) <- output$summary[nn,1:2]
    if (getEdges && is.null(ncores)){
      output$edges[[nn]] <- getFactoredEdges(x,phca$tree)
    }
    
  }
  
  
  if (getEdges && !is.null(ncores)){
    output$edges <- getFactoredEdgesPAR(ncores=ncores,tree=phca$tree,V=phca$basis[,1:ncomponents,drop=F])
  }
  

  print(output$summary)
  return(output)
}