#' Internal function for prediction of data matrix given a method of amalgamation and its expectation
#' @export
#' @param yhat Estimated coordinate (expected value of coordinate) from regression on the amalgamation of Grp
#' @param Grp Two element list containing a group and its complement for amalgamation & contrast
#' @param n number of parts in entire dataset
#' @param method Method for amalgamation.
#' @param Pbasis coming soon

PredictAmalgam <- function(yhat,Grp,n){
  #predicts Data matrix given yhat, Grp and amalgamation method
    v <- ilrvec(Grp,n) %>% c(rep(0,n)) %>% matrix(ncol=2,byrow=F)
    prediction <- yhat %>% c(rep(0,length(yhat))) %>% matrix(ncol=2,byrow=F) %>% compositions::ilrInv(v) %>% t %>% matrix(ncol=length(yhat))

  return(prediction)
}
