#' Internal function for prediction of data matrix given a method of amalgamation and its expectation
#' @export
#' @param yhat Estimated coordinate (expected value of coordinate) from regression on the amalgamation of Grp
#' @param Grp Two element list containing a group and its complement for amalgamation & contrast
#' @param n number of parts in entire dataset
#' @param method Method for amalgamation.
#' @param Pbasis coming soon

PredictAmalgam <- function(yhat,Grp,n,method,Pbasis=1){
  #predicts Data matrix given yhat, Grp and amalgamation method
  if (method=='add'){
    ### In this case, our yhats corresopnd to log-ratios of groups' relative abundances, log(grp1/grp2)
    if (is.null(Pbasis)){#our Grp contains all taxa in our set, 1:n.

    } else { #our group does not contain all taxa in set 1:n and we need

    }
    # prediction <- d
  } else { #use ILR approach
    v <- ilrvec(Grp,n) %>% c(rep(0,n)) %>% matrix(ncol=2,byrow=F)
    prediction <- yhat %>% c(rep(0,length(yhat))) %>% matrix(ncol=2,byrow=F) %>% ilrInv(v) %>% t %>% matrix(ncol=length(yhat))
  }

  return(prediction)
}
