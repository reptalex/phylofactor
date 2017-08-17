#' Concise summarization of phylofactor summary
#'
#' @export
#' @param smry PF summary class object. see \code{\link{pf.summary}}
#' @return list containing concise taxa IDs of groups split, as well as their predicted ratios of abundances in different conditions
#' @examples
#' library(phylofactor)
#' data('FTmicrobiome')
#'
#' smry <- pf.summary(FTmicrobiome$PF,FTmicrobiome$taxonomy,factor=1)
#'
#' pf.tidy(smry)

pf.tidy <- function(smry){
    
    ### Be concise about taxa split
    l <- list(unique(smry$TaxaSplit[[1]][,2]),unique(smry$TaxaSplit[[2]][,2]))
    Monophy <- c(smry$group1$is.monophyletic,smry$group2$is.monophyletic)+1
    names(l) <- mapply(list('group1','group2'),FUN=function(y,ix) paste(y,c(', Paraphyletic',', Monophyletic')[ix],sep=''),ix=Monophy)
    i=2

    ### Be concise about their predicted effects ##
    
    if (!is.null(smry$glm)){
      coef <- smry$glm$coefficients
      l <- c(l,list(coef))
      i=i+1
      names(l)[i] <- 'Coefficients'
      
      pp <- predict(smry$glm)
    
      r <- dim(smry$group1$IDs)[1]
      s <- dim(smry$group2$IDs)[1]
      ratios <- exp(sqrt((r+s)/(r*s))*pp)
      
      
  
      l <- c(l,list(ratios))
      i=i+1
      names(l)[i] <- 'Predicted ratio of group1/group2'
    }

    obsRatios <- gMean(smry$group1$otuData)/gMean(smry$group2$otuData)
    if (!is.null(smry$glm)){
      names(obsRatios) <- names(smry$glm$linear.predictors)
    } else {
      names(obsRatios) <- colnames(smry$group1$otuData)
    }

    l <- c(l,list(obsRatios))
    i=i+1
    names(l)[i] <- 'Observed Ratio of group1/group2 geometric means'

    return(l)
}
