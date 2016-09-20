#' amalgamate two groups in a compositional data matrix to obtain ILR vector
#' Must input the log of the data matrix
#' @export
#' @param Grp List containing two non-overlapping vectors of positive, non-zero integers less than the number of rows in the Data matrix.
#' @param Data Compositional data matrix whose rows are parts and whose columns are samples.

amalg.ILR <- function(Grp,LogData){
    r = length(Grp[[1]])
    s = length(Grp[[2]])
    if (r>1){
      output <- colSums(LogData[Grp[[1]],])*(sqrt(s/(r*(r+s))))
    } else {
      output <- LogData[Grp[[1]],]*(sqrt(s/(r*(r+s))))
    }
    if (s>1){
      output <- output-colSums(LogData[Grp[[2]],])*sqrt(r/(s*(r+s)))
    } else {
      output <- output-LogData[Grp[[2]],]*sqrt(r/(s*(r+s)))
    }
    
    return(output)
}