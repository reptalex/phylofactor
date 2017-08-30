#' amalgamate two groups in a compositional data matrix to obtain ILR vector
#' Must input the log of the data matrix
#' @export
#' @param Grp List containing two non-overlapping vectors of positive, non-zero integers less than the number of rows in the Data matrix.
#' @param LogData Logarithm of data matrix whose rows are parts and whose columns are samples.

amalg.ILR <- function(Grp,LogData){
    # r = length(Grp[[1]])
    # s = length(Grp[[2]])
    # # # if (r>1){
    # # #  output <- colSums(LogData[Grp[[1]],,drop=F])*(sqrt(s/(r*(r+s))))
    # # # } else {
    # # #  output <- LogData[Grp[[1]],,drop=F]*(sqrt(s/(r*(r+s))))
    # # # }
    # # # if (s>1){
    # # #  output <- output-colSums(LogData[Grp[[2]],,drop=F])*sqrt(r/(s*(r+s)))
    # # # } else {
    # # #  output <- output-LogData[Grp[[2]],,drop=F]*sqrt(r/(s*(r+s)))
    # # # }
    # # ############ it's slightly more numerically stable to use the computation below:
    # myColMean <- function(A){
    #   myMean <- function(x) {
    #     mn=0
    #     for (i in 1:length(x)){
    #       mn=mn+(1/i)*(x[i]-mn)
    #     }
    #     return(mn)
    #   }
    #   
    #   return(apply(A,2,FUN=myMean))
    # }
    r = length(Grp[[1]])
    s = length(Grp[[2]])
    output <- (colMeans(LogData[Grp[[1]],,drop=F])-colMeans(LogData[Grp[[2]],,drop=F]))*sqrt(r*s/(r+s))
    # output <- (myColMean(LogData[Grp[[1]],,drop=F])-myColMean(LogData[Grp[[2]],,drop=F]))*sqrt(r*s/(r+s))
    ############# the most numerically stable option is in the compositions package
    ############# but for some reason the class of object returned screws up the rest of phylofactor
    # LogData <- compositions::rcomp(t(Data))
    # v <- ilrvec(Grp,ncol(LogData)) %>% matrix(.,ncol=1)
    # output <- suppressWarnings(apply(LogData,2,FUN=function(x,v) compositions::ilr(,v),v=v))
    return(output)
}