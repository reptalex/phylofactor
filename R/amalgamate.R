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
      rownames(output) <- mapply(FUN = paste, as.list(repa('Atom',n)),as.list(1:n))
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
