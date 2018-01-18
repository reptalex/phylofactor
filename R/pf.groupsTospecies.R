#' Converts integer-indexed groups to species-labelled groups in phylofactor object
#' 
#' @export
#' @param pf phylofactor class object from \code{\link{PhyloFactor}}, \code{\link{twoSampleFactor}} or \code{\link{gpf}}
#' @examples 
#' library(phylofactor)
#' data(FTmicrobiome)
#' pf <- FTmicrobiome$PF
#' species.groups <- pf.groupsTospecies(pf)
#' pf$factors[1:3,]
#' #the otuIDs/species in factor 3, Group1 are:
#' species.groups[[3]][[1]]
pf.groupsTospecies <- function(pf){
  return(lapply(pf$groups,FUN=function(grp,spp) lapply(grp,FUN=function(grp,spp) spp[grp],spp=spp),spp=pf$tree$tip.label))
}