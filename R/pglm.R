#' Internal wrapper for \code{\link{glm}} for phylofactor
#' @export
#' @param y response variable
#' @param xx independent variable
#' @param frmla Formula for dependence of y on x
#' @param smallglm Logical. See \code{\link{PhyloFactor}}
#' @param ... optional input arguments to \code{\link{glm}}
#' @return glm object
#' @examples
#' x <- rnorm(10)
#' y <- x+rnorm(10)
#' pglm(y,x,y~x)

################# phyloreg ############################
pglm <- function(y,xx,frmla,smallglm=T,...){
  ##performs individual regression.
  dataset <- c(list(y),as.list(xx))
  names(dataset) <- c('Data',names(xx))
  dataset <- model.frame(frmla,data = dataset)
  
  if(smallglm){
    return(glm(frmla,data = dataset,...))
  } else {
    return(biglm::bigglm(frmla,data = dataset,...))
  }
}
