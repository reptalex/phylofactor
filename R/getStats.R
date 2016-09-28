#' Grabs P value and F statistic from aov of glm object.
#' @export
#' @param GLM glm object
#' @param y dependent variable - used when inputting bigglm object.
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
  if (!(class(GLM)[1] %in% c('lm','gam','glm','bigglm'))){stop('unknown input class - need glm or bigglm')}

  if (class(GLM)[1] %in% c('lm','gam','glm')){
    av <- aov(GLM)
    P <- summary(av)[[1]][1,5]
    f <- summary(av)[[1]][1,4]
    ex.var <- (GLM$null.deviance-GLM$deviance)/length(GLM$residuals)
    stats <- c(P,f,ex.var)
    names(stats) <- c('Pval','F','ExplainedVar')
  } else {
    if (is.null(y)){stop('if inputting bigglm, need to input dependent variable, y')}
    dfr <- GLM$df.resid
    dfx <- (GLM$n-1-dfr)
    ssq <- sum((y-mean(y))^2)-GLM$rss
    msq <- ssq/dfx
    rss <- GLM$rss
    mrss <- rss/GLM$df.resid
    
    ex.var <- ssq/(GLM$n-1)

    Fstat <- msq/mrss
    Pval <-  1-pf(Fstat,df1=dfx,df2=dfr)

    stats <- c(Pval,Fstat,ex.var)
    names(stats) <- c('Pval','F','ExplainedVar')
  }
  return(stats)
}
