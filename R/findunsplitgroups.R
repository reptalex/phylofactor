#' Finds unsplit groups in a signed partition matrix, V
#'
#' @param V Partition matrix.
#' @example
#' V <- matrix(c(1,1,-1,-1,1,-1,0,0),byrow=F,ncol=2)
#' find.unsplit.Grps(V)

#################################### find.unsplit.groups ###############################################
#useful function for grabbing groups given a partition

find.unsplit.Grps <- function(V){ #function to find out the groups which still need to be partitioned
  #input ilr sub-basis, V. output a list of unsplit groups.
  S <- sign(V)
  if (sum(S[,1]==0)>0){stop('obtaining atoms for sub-trees without the initial partition are currently unavailable. Try inputting the whole basis up to your desired clade')}
  if (is.null(dim(S))){
    ix <- 1:length(S)
    m <- 1
    ss <- vector(mode='list',length=m+1)
    ss[1:2] <- split(ix,S)
  } else {
    ix=1:dim(S)[1]
    m <- dim(S)[2]
    ss <- vector(mode='list',length=m+1)
    ss[1:length(unique(S[,1]))] <- split(ix,S[,1])
  }


  if(m>1){
    for (ll in 2:m){
      ind <- which(S[,ll]!=0)
      inx <- which(unlist(lapply(ss,function(x,y){all(y %in% x && x %in% y)},y=ind)))
      dum <- split(ss[[inx]],S[which(S[,ll]!=0),ll])
      ss[[inx[1]]]=dum[[1]]
      ss[[ll+1]]=dum[[2]]
    }}
  ss <- ss[which(unlist(lapply(ss,function(x){return(length(x)>1)})))]
  return(ss)
}
