#' Internal wrapper for \code{\link{glm}} for phylofactor
#' @export
#' @param y response variable
#' @param xx independent variable
#' @param frmla Formula for dependence of y on x
#' @param ... optional input arguments to \code{\link{glm}}
#' @return glm object
#' @examples
#' x <- rnorm(10)
#' y <- x+rnorm(10)
#' pglm(y,x,y~x)

################# phyloreg ############################
pglm <- function(Y,xx,frmla,...){
  ##performs individual regression.
  dataset <- c(list(Y),as.list(xx))
  names(dataset) <- c('Data',names(xx))
  args <- list('data'=dataset,'formula'=frmla,...)
  return(do.call(stats::glm,args))
}
