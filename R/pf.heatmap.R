#' Heatmap for phylofactor data
#' @export
#' @param PF phylofactor class object from either \code{\link{PhyloFactor}} or \code{\link{gpf}} with \code{algorithm='mStable'}
#' @param tree optional input, only used if \code{PF} is not input. If inputting \code{tree}, must also input \code{Data}
#' @param Data Dataset with rows equal to the number of species. Ideally, rownames are all equal to \code{PF$tree} of \code{tree} tip-labels. Otherwise, rows must be in order of tree tip-labels.
#' @param factors vector of integer factors for input with \code{PF}. Will use \code{\link{pf.tree}} to highlight clades obtained at each factor.
#' @param column.order vector of integers of length equal to \code{ncol(Data)} or \code{ncol(PF$Data)}. Will re-order the data (e.g. according to desired meta-data).
#' @param branch.length input for \code{\link{ggtree}}
#' @param color input \code{color} for \code{gheatmap}
#' @param font.size font size for \code{gheatmap}
#' @param ... additional arguments passed to \code{\link{gheatmap}}
#' @examples
#' library(phylofactor)
#' require(ggpubr)
#' data(FTmicrobiome)
#' PF <- FTmicrobiome$PF
#' obs <- PF$transform.fcn(PF$Data)
#' pred <- predict(PF,factor=3)
#' observed <- pf.heatmap(PF,Data=obs,factors=1:3,width=3)
#' predicted <- pf.heatmap(PF,Data=pred,factors=1:3,width=3)
#' ggarrange(observed,predicted,nrow=2)

pf.heatmap <- function(PF=NULL,tree=NULL,Data=NULL,factors=NULL,column.order=NULL,branch.length='none',color=NA,font.size=0,...){
  
  if (is.null(PF)){
    if (is.null(tree) | is.null(Data)){
      stop('If not inputting phylofactor object, must input tree & data')
    } else {
      PF <- NULL
      PF$tree <- tree
    }
    if (!is.null(factors)){
      stop('cannot input factors without input PF phylofactor object')
    }
  } else {
    if (!PF$phylofactor.fcn %in% c('PhyloFactor','gpf')){
      stop('pf.heatmap only works for PhyloFactor or gpf mStable input')
    }
    if (PF$phylofactor.fcn=='gpf'){
      if (PF$algorithm!='mStable'){
        stop('gpf-based pf.heatmap only works for algorithm==mStable')
      } else {
        if (class(PF$Data)=='list'){
          if (is.null(Data)){
            Data <- PF$Data$Successes/(PF$Data$Successes+PF$Data$Failures)
          }
        }
      }
    }
  }
  
  
  if (!is.null(factors)){
    if (max(factors)>PF$nfactors){
      stop('max of factors cannot exceed PF$nfactors')
    }
  }
  
  if (!is.null(Data)){
    if (nrow(Data) != length(PF$tree$tip.label)){
      stop('number of rows of Data does not equal number of tips in PF$tree')
    } else {
      if (is.null(rownames(Data))){
        warning('No rownames of input Data - pf.heatmap will assume rows are in same order of tree tip labels')
        rownames(Data) <- PF$tree$tip.label
      } else {
        if (!all(rownames(Data) %in% PF$tree$tip.label & all(PF$tree$tip.label %in% rownames(Data)))){
            stop('could not match all rownames of Data with all tree tip labels.')
        } else {
            Data <- Data[PF$tree$tip.label,]
        }
      }
    }
  }
  
  if (is.null(Data)){
      Data <- as.data.frame(PF$transform.fcn(PF$Data))
  } else {
    if (class(Data)!='data.frame'){
      nms <- colnames(Data)
      Data <- as.data.frame(Data)
      if (is.null(nms)){
        colnames(Data) <- 1:ncol(Data)
      }
    }
  }
  
  if (is.null(factors)){
    gg <- ggtree::ggtree(PF$tree,layout='rectangular',branch.length = branch.length)
  } else {
    suppressWarnings(gg <- pf.tree(PF,layout='rectangular',factors = factors,branch.length=branch.length)$ggplot)
  }
  if (is.null(column.order)){
    column.order <- 1:ncol(Data)
  }
  gg <- ggtree::gheatmap(gg,Data[,column.order],color=color,font.size=font.size,...)
  return(gg)
}
