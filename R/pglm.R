#' Returns regression of x on y according to frmla

#' @param y response variable
#' @param x independent variable
#' @param frmla Formula for dependence of y on x
#' @return glm object
#' @examples
#' x <- rnorm(10)
#' y <- x+rnorm(10)
#' phyloreg(y,x,y~x)

################# phyloreg ############################
pglm <- function(y,x,frmla,...){
  #input: dependent variable y, independent variable x, formula 'frmla' and 'choice' to be used to
  #compare different phylogenetic partitions
#   ##performs individual regression.
#   M <- data.frame(y,x)
#   names(M) <- c('Data','X')
  return(glm(frmla,data = data.frame('Data'=y,'X'=x),...))
}
