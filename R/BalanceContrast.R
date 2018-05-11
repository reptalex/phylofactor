#' Project data matrix onto balance contrast (if data are log-data, then this will produce ILR balance)
#' @export
#' @param Grp List containing two non-overlapping vectors of positive, non-zero integers less than the number of rows in the Data matrix.
#' @param TransformedData Transformed data matrix whose rows are parts and whose columns are samples.

BalanceContrast <- function(Grp,TransformedData){
   
    r = length(Grp[[1]])
    s = length(Grp[[2]])
    output <- (colMeans(TransformedData[Grp[[1]],,drop=F])-colMeans(TransformedData[Grp[[2]],,drop=F]))*sqrt(r*s/(r+s))
   
    return(output)
}