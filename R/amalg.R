#' amalgamate two groups in a compositional data matrix to obtain ILR vector
#' Must input the log of the data matrix
#' @export
#' @param Grp List containing two non-overlapping vectors of positive, non-zero integers less than the number of rows in the Data matrix.
#' @param LogData Logarithm of data matrix whose rows are parts and whose columns are samples.

amalg.ILR <- function(Grp,LogData){
   
    r = length(Grp[[1]])
    s = length(Grp[[2]])
    output <- (colMeans(LogData[Grp[[1]],,drop=F])-colMeans(LogData[Grp[[2]],,drop=F]))*sqrt(r*s/(r+s))
   
    return(output)
}