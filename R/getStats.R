#' Grabs P value, F statistic and difference between null and residual deviance from glm, lm, gam and bigglm objects.
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
  if (!(any(c('lm','gam','glm','bigglm') %in% class(GLM)))){stop('unknown input class - need glm or bigglm')}

  if (!'bigglm' %in% class(GLM)){
    dfr <- GLM$df.residual
    dfx <- GLM$df.null-dfr
    ssq <- (GLM$null.deviance-GLM$deviance)
    msq <- ssq/dfx
    mrss <- GLM$deviance/GLM$df.residual
    
    f <- msq/mrss
    P <- 1-pf(f,df1=dfx,df2=dfr)
    ex.var <- ssq/length(GLM$residuals)
    stats <- c(P,f,ex.var)
    
    names(stats) <- c('Pval','F','ExplainedVar')
  } else {
    if (is.null(y)){stop('if inputting bigglm, need to input dependent variable, y')}
    dfr <- GLM$df.resid
    dfx <- (GLM$n-1-dfr)
    ssx <- sum((y-mean(y))^2)-GLM$rss
    msq <- ssx/dfx
    mrss <- GLM$rss/GLM$df.resid
    
    ex.var <- ssx/(GLM$n-1)

    Fstat <- msq/mrss
    Pval <-  1-pf(Fstat,df1=dfx,df2=dfr)

    stats <- c(Pval,Fstat,ex.var)
    names(stats) <- c('Pval','F','ExplainedVar')
  }
  return(stats)
}
