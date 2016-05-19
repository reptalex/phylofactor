#' Summary of phylofactor object for a given node or factor number
#' @export
#' @param PF PhyloFactor object
#' @param tree phylogeny (ape class)
#' @param taxonomy Taxonomy, first column is OTU ids in tree, second column is greengenes taxonomic string
#' @param node node number, must be in PF$nodes. Takes priority over "factor" for what is summarized
#' @param factor Factor number to summarize.
#' @param simplify.TaxaIDs Whether or not to simplify IDs of group taxonomy to their shorest unique prefix (for each OTU in each group, this is taxonomy down to the coarsest taxonomic level unique to the OTU's group)
#' @param subtree Logical indicating whether or not to output the subtree partitioned when PF$node was split (i.e. sub-tree split by factor)
#' @param prediction Logical. If subtree=T, prediction=T will produce a phylo.heatmap containing the predicted data. Otherwise, will output phylo.heatmap of real data
#' @param tipLabels Logical indicating whether or not to include tip labels in plot of subtree.
#' @param Transform String, method for transforming data for easier visualization of sub-tree phylo.heatmaps. Options include "atan" to flatten extreme values, "tan" to amplify extreme values, and "clr" for centered log-ratio transform. 
#' @param ... additional arguments for phylo.Heatmap. If phyloHeatmap='AW', these are additional arguments for \code{plot.phylo}. Otherwise, see description for \code{phylo.heatmap} in phytools package for list of additional arguments.
#' @return summary object. List containing $group and $complement info, each containing summary.group output for that group -  $IDs, $otuData and $PF.prediction
#' @examples
#' data("FTmicrobiome")
#' OTUTable <- FTmicrobiome$OTUTable        #OTU table
#' Taxonomy <- FTmicrobiome$taxonomy        #taxonomy
#' tree <- FTmicrobiome$tree                #tree
#' X <- FTmicrobiome$X                      #independent variable - factor indicating if sample is from feces or tongue
#'
#' rm('FTmicrobiome')
#'
#' # remove rare taxa
#' ix <- which(rowSums(OTUTable==0)<30)
#' OTUTable <- OTUTable[ix,]
#' OTUs <- rownames(OTUTable)
#' tree <- drop.tip(tree,which(!(tree$tip.label %in% OTUs)))
#'
#' par(mfrow=c(1,1))
#' phylo.heatmap(tree,t(clr(t(OTUTable))))
#' PF <- PhyloFactor(OTUTable,tree,X,nfactors=2,choice='var')
#'
#' FactorSummary <- phylofactor.summary(PF,Taxonomy,factor=1)
#'
#' str(FactorSummary)
#'
#' par(mfrow=c(1,2))
#' plot(FactorSummary$ilr,ylab='ILR coordinate',main='ILR coordinate of factor',xlab='sample no.',pch=16)
#' lines(FactorSummary$fitted.values,lwd=2,col='blue')
#' legend(x=1,y=-5,list('data','prediction'),pch=c(16,NA),lty=c(NA,1),col=c('black','blue'),lwd=c(NA,2))
#'
#' plot(FactorSummary$MeanRatio,ylab='ILR coordinate',main='Mean Ratio of Grp1/Grp2',xlab='sample no.',pch=16)
#' lines(FactorSummary$fittedMeanRatio,lwd=2,col='blue')
#' legend(x=1,y=-5,list('data','prediction'),pch=c(16,NA),lty=c(NA,1),col=c('black','blue'),lwd=c(NA,2))
phylofactor.summary <- function(PF,taxonomy,factor=NULL,tree=PF$tree,simplify.TaxaIDs=F,subtree=F,prediction=T,tipLabels=F,Transform="atan",...){
  #summarizes the IDs of taxa for a given node identified as important by PhyloFactor. If subtree==T, will also plot a subtree showing the taxa
  if (is.null(factor)){stop('need to input a factor')}
  grp1 <- PF$groups[[factor]][[1]]
  grp2 <- PF$groups[[factor]][[2]]

  output <- NULL
  output$group1 <- summary.group(PF,tree,taxonomy,factor = factor,grp1,simplify=simplify.TaxaIDs)
  output$group2 <- summary.group(PF,tree,taxonomy,factor = factor,grp2,simplify=simplify.TaxaIDs)



  if (subtree==T){
    if (factor>1){
      atms <- atoms(PF$basis[,1:(factor-1),drop=F])
      grp <- c(grp1,grp2)
      tree <- ape::drop.tip(tree,setdiff(unlist(atms),grp))
    }
    edgs <- rep(2,ape::Nedge(tree))
    cols <- rep('black',ape::Nedge(tree))

    edgG <- extractEdges(tree,taxa=tree$tip.label[grp1],type=3)
    edgs[edgG]=8
    cols[edgG]='red'
      if (prediction==T){
        dta <- rbind(output$group1$PF.prediction,output$group2$PF.prediction)
      } else { 
        dta <- t(compositions::clr(t(rbind(output$group1$otuData,output$group2$otuData))))
      }
    if (Transform=='atan'){
      dta <- atan(dta)
    } else {
      dta <- t(compositions::clr(t(dta)))
    }
    phylo.heatmapAW(tree=tree,Y=dta,tipLabels = tipLabels,edge.width=edgs,edge.color=cols,...)
    output$subtree <- tree
  }

  output$TaxaSplit <- TaxaSplit(output)
  output$glm <- PF$glms[[factor]]
  output$ilr <- PF$glms[[factor]]$y
  output$fitted.values <- PF$glms[[factor]]$fitted.values

  r <- length(grp1)
  s <- length(grp2)
  output$MeanRatio <- exp(output$ilr/(sqrt(r*s/(r+s))))
  output$fittedMeanRatio <- exp(output$fitted.values/(sqrt(r*s/(r+s))))


  class(output) <- 'PF summary'
  return(output)

}

