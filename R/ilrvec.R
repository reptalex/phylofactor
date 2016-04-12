#' Internal function for computing an ilr balancing element given a list of split taxa in and the total number of taxa
#'
#' @param ix two-member list containing the parts split by a partition
#' @param n total number of parts
#' @return ILR basis element corresponding to the split of ix[[1]] with ix[[2]] out of n-part composition.
#' @examples
#' n <- 10
#' ix <- NULL
#' ix[[1]] <- c(3,4)
#' ix[[2]] <- c(5,6,7,8)
#' v1 <- ilrvec(ix,n)
#' sum(v1)
#'
#' ix[[1]] <- c(5,6)
#' ix[[2]] <- c(7,8)
#' v2 <- ilrvec(ix,n)
#' sum(v2)
#' v1 %*% v2

######################### ilrvec ###############################
ilrvec <- function(ix,n){
  # Computes ilr basis element for n-species composition corresponding to group given by 2-element list, "ix".
  # Input indexes, ix, which are partitioned.

  if (!all(unlist(ix) %in% 1:n)){stop('indexes in ix are not in 1:n')}
  if (any(ix[[1]] %in% ix[[2]])){stop('indexes ix have non-empty intersect')}

    ix1=unlist(ix[[1]])
    ix2=unlist(ix[[2]])
    r <- length(ix1)
    s <- length(ix2)
    a <- sqrt(s/(r*(r+s))) #ix1 elements
    b <- -sqrt(r/(s*(r+s))) #ix2 elements
    vec=rep(0,n)
    vec[ix1]=a
    vec[ix2]=b

  return(vec)
}
