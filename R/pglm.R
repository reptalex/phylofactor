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
  #input: dependent variable y, independent variable x, formula 'frmla' and 'choice' to be used to
  #compare different phylogenetic partitions
#   ##performs individual regression.
#   M <- data.frame(y,x)
#   names(M) <- c('Data','X')
  dataset <- model.frame(frmla,data = list('Data'=y,'X'=xx))
  if(smallglm){
    return(glm(frmla,data = dataset,...))
  } else {
    return(biglm::bigglm(frmla,data = dataset,...))
  }
}
