#' amalgamate two groups in a compositional data matrix
#' @export
#' @param Grp List containing two non-overlapping vectors of positive, non-zero integers less than the number of rows in the Data matrix.
#' @param Data Compositional data matrix whose rows are parts and whose columns are samples.
#' @param method Method for amalgamation. Takes values of "sum" and "ILR". ILR amalgamates using geometric means, whereas "sum" amalgamates by summing relative abundances.
#' @param collapse Logical indicating whether or not to collapse groups (collapse=T) by summing or geometric means, or to output log-ratios from amalgamation & comparison.
#' @return If collapse=T, returns either a log-ratio of the groups. If collapse=F, returns the geometric mean of each group.
#' @examples
#' library(compositions)
#' set.seed(1)
#' Data <- matrix(rlnorm(80),nrow=8)
#' Grp <- list(c(1,2),c(3,4,8))
#' for (gg in Grp[[1]]){
#' Data[gg,] <- Data[gg,]*10^(seq(-1,1,length.out=10))
#' }
#'
#' Data <- t(rcomp(t(Data)))
#'
#' method='ILR'
#' lr <- amalgamate(Grp,Data,method)
#' plot(1:10,lr,xlab='sample',ylab='log-ratio of Grp {1,2} over {3,4,8}')

amalgamate <- function(Grp,Data,method,collapse=F){


  if (collapse==T){

      n <- length(Grp) #this is the number of rows in our new OTU table
      m <- dim(OTUTable)[2]
      output <- matrix(0,n,m)

      for (nn in 1:n){
        grp <- Grp[[nn]]
          if (method=='add'){
            if (length(grp)>1){
              output[nn,] <- geometricmeanCol(Data[grp,])
            } else {
              output[nn,] <- Data[grp,]
            }
          } else {
            if (length(grp)>1){
              output[nn,] <- colSums(Data[grp,])
            } else {
              output[nn,] <- Data[grp,]
            }
          }
      }

      output <- output %>% t %>% clo %>% t
      rownames(output) <- mapply(FUN = paste, as.list(rep('Atom',n)),as.list(1:n))
      return(output)

    } else {

    grp1=Grp[[1]]
    r = length(grp1)
    grp2=Grp[[2]]
    s = length(grp2)

    if (method=='add'){
      if (r>1 && s==1){
        return(log(colSums(Data[grp1,])/Data[grp2,]))
      }
      if (r==1 && s>1){
        return(log(Data[grp1,]/colSums(Data[grp2,])))
      }
      if (r==1 && s==1){
        return(log(Data[grp1,]/Data[grp2,]))
      }
      if (r > 1 && s>1){
        return(log(colSums(Data[grp1,])/colSums(Data[grp2,])))
      }
    } else {
      if (r>1){
        b1 <- apply(Data[grp1,],MARGIN = 2,FUN  = function(x){return(sum(log(x)))})*(sqrt(s/(r*(r+s))))
      } else {
        b1 <- log(Data[grp1,])*sqrt(s/(r*(r+s)))
      }
      if (s>1){
        b2 <- apply(Data[grp2,],MARGIN = 2,FUN = function(x){return(sum(log(x)))})*sqrt(r/(s*(r+s)))
      } else {
        b2 <- log(Data[grp2,])*sqrt(r/(s*(r+s)))
      }
      return(b1-b2)
    }
  }
}
