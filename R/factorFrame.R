#' Function for making factor-contrast in \code{\link{gpf}}
#' @export
#' @param grp list containing two disjoint lists of species, such as thouse output from \code{\link{getGroups}}
#' @param spp set of species or tip-labels of tree.
factorFrame <- function(grp,spp){
  g1 <- spp[grp[[1]]]
  r <- length(g1)
  g2 <- spp[grp[[2]]]
  s <- length(g2)
  nas <- setdiff(spp,spp[unlist(grp)])
  n<- length(nas)
  ff <- data.table('Species'=c(g1,g2,nas),'Size'=c(rep(r,r),rep(s,s),rep(NA,n)),'G'=factor(c(rep('R',r),rep('S',s),rep(NA,n))),key='Species')
  rownames(ff) <- ff$Species
  return(ff)
}
