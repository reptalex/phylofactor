#' internal function for phylofactor
preg <- function(y,xx,frmla,RegressionFunction,...){
  dataset=as.data.frame(cbind(y,xx))
  if (is.null(dim(X))){
    names(dataset) <- c('Data','X')
  } else {
    names(dataset) <- c('Data',colnames(xx))
  }
  dataset <- model.frame(frmla,data=dataset)
  return(RegressionFunction(frmla,data=dataset,...))
}