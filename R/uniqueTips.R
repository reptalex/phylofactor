uniqueTips <- function(tree,sigClades){
  Clds <- length(sigClades)
  tips <- NULL
  for (cc in 1:Clds){
    tips <- sort(unique(c(tips,Descendants(tree,sigClades[cc],type='tips')[[1]])))
  }
  return(tips)
}
