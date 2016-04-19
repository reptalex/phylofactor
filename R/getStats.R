#' Grabs P value and F statistic from aov of glm object.
#' @export
#' @param GLM glm object
#' @return vector containing P value and F statistic from aov(GLM)
#' @examples
#' x <- rnorm(10)
#' y <- 2*x+rnorm(10)
#'
#' GLM <- glm(y~x)
#' getStats(GLM)

############################# getP #####################################
getStats <- function(GLM,y=NULL){
  #extract P value and F statistic from analysis of variance for an input glm
  if (!(class(GLM)[1] %in% c('glm','bigglm'))){stop('unknown input class - need glm or bigglm')}

  if (class(GLM)[1]=='glm'){
    av <- aov(GLM)
    P <- summary(av)[[1]][1,5]
    f <- summary(av)[[1]][1,4]
    stats <- c(P,f)
    names(stats) <- c('Pval','F')
  } else {
    if (is.null(y)){stop('if inputting bigglm, need to input data, y')}
    dfx <- (GLM$n-GLM$df.resid)
    dfr <- GLM$df.resid
    ssq <- sum(y^2)-GLM$rss
    msq <- ssq/dfx
    rss <- GLM$rss
    mrss <- rss/GLM$df.resid

    Fstat <- msq/mrss
    Pval <-  1-pf(Fstat,df1=dfx,df2=dfr)

    stats <- c(Pval,Fstat)
    names(stats) <- c('Pval','F')
  }
  return(stats)
}
