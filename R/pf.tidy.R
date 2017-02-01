#' Concise summarization of phylofactor summary
#'
#' @export
#' @param smry PF summary class object. see \code{\link{phylofactor.summary}}
#' @return list containing concise taxa IDs of groups split, as well as their predicted ratios of abundances in different conditions
#' @examples
#' library(phylofactor)
#' data('FTmicrobiome')
#'
#' smry <- phylofactor.summary(FTmicrobiome$PF,FTmicrobiome$taxonomy,factor=1)
#'
#' pf.tidy(smry)

pf.tidy <- function(smry){
    
    ### Be concise about taxa split
    l <- list(unique(smry$TaxaSplit[[1]][,2]),unique(smry$TaxaSplit[[2]][,2]))
    Monophy <- c(smry$group1$is.monophyletic,smry$group2$is.monophyletic)+1
    names(l) <- mapply(list('group1','group2'),FUN=function(y,ix) paste(y,c(', Paraphyletic',', Monophyletic')[ix],sep=''),ix=Monophy)


    ### Be concise about their predicted effects ##
    ### If 
      coef <- smry$glm$coefficients
      l <- c(l,list(coef))
      names(l)[3] <- 'Coefficients'
      
      pp <- predict(smry$glm)
    
    r <- dim(smry$group1$IDs)[1]
    s <- dim(smry$group2$IDs)[1]
    ratios <- exp(sqrt((r+s)/(r*s))*pp)
    
    

    l <- c(l,list(ratios))
    names(l)[4] <- 'Predicted ratio of group1/group2'


    obsRatios <- gMean(smry$group1$otuData)/gMean(smry$group2$otuData)
    names(obsRatios) <- names(smry$glm$linear.predictors)

    l <- c(l,list(obsRatios))
    names(l)[5] <- 'Observed Ratio of group1/group2 geometric means'

    return(l)
}
