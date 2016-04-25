#' Internal function - id which atom in atomList split
#' @param grp split atom
#' @param atomList list of atoms
#' @return index where atomList[[ix]] contains grp
whichAtomSplit <- function(grp,atomList){
  return(which(sapply(atomList,FUN=function(atm,g) all(g %in% atm) && all(atm %in% g),g=unlist(grp),simplify=T)))
}
