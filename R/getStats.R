############################# getP #####################################
getStats <- function(gg){
  #extract P value and F statistic from analysis of variance for an input glm
  av <- aov(gg)
  P <- summary(av)[[1]][1,5]
  f <- summary(av)[[1]][1,4]
  stats <- c(P,f)
  names(stats) <- c('Pval','F')
  return(stats)
}