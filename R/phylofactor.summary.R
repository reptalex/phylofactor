#' Summary of phylofactor object for a given node or factor number
#' @export
#' @param PF PhyloFactor object
#' @param tree phylogeny (ape class)
#' @param taxonomy Taxonomy, first column is OTU ids in tree, second column is greengenes taxonomic string
#' @param node node number, must be in PF$nodes. Takes priority over "factor" for what is summarized
#' @param factor Factor number to summarize.
#' @param subtree Logical indicating whether or not to output the subtree partitioned when PF$node was split (i.e. sub-tree split by factor)
#' @param prediction Logical. If subtree=T, prediction=T will produce a phylo.heatmap containing the predicted data. Otherwise, will output phylo.heatmap of real data
#' @param tipLabels Logical indicating whether or not to include tip labels in plot of subtree.
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
#' delta=0.65
#' OTUTable[OTUTable==0]=delta
#' OTUTable <- OTUTable %>% t %>% clo %>% t
#' OTUTable <- OTUTable[tree$tip.label,]
#'
#' par(mfrow=c(1,1))
#' phylo.heatmap(tree,t(clr(t(OTUTable))))
#' PF <- PhyloFactor(OTUTable,tree,X,nclades=2,choice='var')
#'
#' FactorSummary <- summary.phylofactor(PF,tree,Taxonomy,factor=1)
#' NodeSummary <- summary.phylofactor(PF,tree,Taxonomy, node=PF$nodes[2])
#'
#' str(FactorSummary)
#'
#' NodeSummary$group$IDs
#'
phylofactor.summary <- function(PF,tree,taxonomy,factor=NULL,subtree=F,prediction=T,tipLabels=F,...){
  #summarizes the IDs of taxa for a given node identified as important by PhyloFactor. If subtree==T, will also plot a subtree showing the taxa
  if (is.null(factor)){stop('need to input a factor')}
  grp1 <- PF$groups[[factor]][[1]]
  grp2 <- PF$groups[[factor]][[2]]

  output <- NULL
  output$group1 <- summary.group(PF,tree,taxonomy,factor = factor,grp1)
  output$group2 <- summary.group(PF,tree,taxonomy,factor=factor,grp2)

  if (subtree==T){
    tr <- ape::drop.tip(tree,setdiff(unlist(atms),grp))
    edgs <- rep(2,ape::Nedge(tr))
    cols <- rep('black',ape::Nedge(tr))

    edgG <- extractEdges(tr,taxa=tree$tip.label[grp1],type=3)
    edgs[edgG]=8
    cols[edgG]='red'
    if (prediction==T){
    dta <- rbind(output$group$PF.prediction,output$complement$PF.prediction)
    } else { dta <- t(compositions::clr(t(rbind(output$group$otuData,output$complement$otuData))))}
    phylo.heatmapAW(tree=tr,Y=dta,tipLabels = tipLabels,edge.width=edgs,edge.color=cols,...)
    output$subtree <- tr
  }

  return(output)

}

