#' Closure function
#' @export
#' @param A matrix or vector
#' @examples
#' clo(c(1,2,3))
clo <- function(A){
  if (is.null(ncol(A))){
    nms <- names(A)
    A <- matrix(A,ncol=length(A))
    rename=T
  } else {
      rename=F
    }
  if (any(A<0)){stop('Compositional matrix input to clo has negative values')}
  output <- A/rowSums(A)
  if (rename){
    colnames(output) <- nms
  }
  return(output)
}