residualVar <- function(prediction,Data){
  return(var(c(t(clr(t(Data))-clr(t(prediction))))))
}
