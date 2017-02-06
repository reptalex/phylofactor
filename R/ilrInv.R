#' Function to take inverse of ilr transform
#' 
#' @export
#' @param x balances
#' @param V balancing elements

ilrInv <- function(x,V){
  x <- t(x)
  M <- t(clo(t((exp(V %*% x)))))
  rownames(M) <- rownames(V)
  if (!is.null(nrow(x))){
    colnames(M) <- colnames(x)
  } else {
    colnames(M) <- names(x)
  }
  return(M)
}