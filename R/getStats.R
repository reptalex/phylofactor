#' Grabs P value and F statistic from aov of glm object.
#'
#' @param GLM glm object
#' @return vector containing P value and F statistic from aov(GLM)
#' @examples
#' x <- rnorm(10)
#' y <- 2*x+rnorm(10)
#'
#' GLM <- glm(y~x)
#' getStats(GLM)

############################# getP #####################################
getStats <- function(GLM){
  #extract P value and F statistic from analysis of variance for an input glm
  av <- aov(GLM)
  P <- summary(av)[[1]][1,5]
  f <- summary(av)[[1]][1,4]
  stats <- c(P,f)
  names(stats) <- c('Pval','F')
  return(stats)
}
