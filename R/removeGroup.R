#' Slow, brute-force AF internal function for removing a group from a list of groups/complements
#' @export
#' @param Grps Full set of groups and their complements, in format of output of \code{\link{getGroups}}
#' @param id Group number to be removed from Grps
#' @return list of groups and their complements constrained to the partition defined by removing Grps[[id]]
#' @examples
#' set.seed(1)
#' tree <- rtree(10)
#' Grps <- getGroups(tree)
#' removeGroup(Grps,11)

removeGroup <- function(Grps,id){
  ##### This code can definitely be made cleaner & more efficient.
  grp <- Grps[[id]]  #find the chosen group
  Grps <- Grps[setdiff(1:length(Grps),id)] #remove chosen group


  #Now put remaining groups into their appropriate partitions
  #instead of contrasts being between G1/G2, they are now (G1\grp)/(G2\grp)
  for (nn in 1:2){
    for (mm in 1:length(Grps)){
      for (kk in 1:2){
        if(all(grp[[nn]] %in% Grps[[mm]][[kk]])){Grps[[mm]][[kk]] <- setdiff(Grps[[mm]][[kk]],grp[[nn]])}
      }
    }
  }

  #some remaining groups will be redundant, namely those groups of clades neighboring the new partition.
  if (length(Grps)>1){
    ix <- NULL
    for (mm in 1:(length(Grps)-1)){
      for(nn in (mm+1):length(Grps)){
        if (all(Grps[[mm]][[1]] %in% Grps[[nn]][[2]])  && all(Grps[[mm]][[2]] %in% Grps[[nn]][[1]])){ix <- c(ix,mm)}
        if (all(Grps[[mm]][[1]] %in% Grps[[nn]][[1]])  && all(Grps[[mm]][[2]] %in% Grps[[nn]][[2]])){ix <- c(ix,mm)}
      }
    }
    Grps <- Grps[setdiff(1:length(Grps),ix)]
  }
  return(Grps)
}
