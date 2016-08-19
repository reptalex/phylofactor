#' Returns regression of x on y according to frmla
#' @export
#' @param y response variable
#' @param x independent variable
#' @param frmla Formula for dependence of y on x
#' @return glm object
#' @examples
#' x <- rnorm(10)
#' y <- x+rnorm(10)
#' phyloreg(y,x,y~x)

################# phyloreg ############################
pglm <- function(y,xx,frmla,smallglm=F,...){
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
