#' Internal function - id which bin in binList split
#' @param grp split bin
#' @param binList list of bins
#' @return index where binList[[ix]] contains grp
whichBinSplit <- function(grp,binList){
  return(which(sapply(binList,FUN=function(atm,g) all(g %in% atm) & all(atm %in% g),g=unlist(grp),simplify=T)))
}
