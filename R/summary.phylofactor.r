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
summary.phylofactor <- function(PF,tree,taxonomy,node=NULL,factor=NULL,subtree=F,prediction=T,tipLabels=F,...){
  #summarizes the IDs of taxa for a given node identified as important by PhyloFactor. If subtree==T, will also plot a subtree showing the taxa
  if (is.null(node)==F){
    if (is.null(factor)==F){stop('need to input either a node or a factor, not both')}
    if ((node %in% PF$nodes)==F){stop('Node is not contained in PhyloFactor')}
  } else {
    if (is.null(factor)){stop('need to input either a node or a factor.')}
    if ((length(PF$nodes))>0){
      if (factor > length(PF$nodes)){stop('factor input exceeds the number of factors in Phylofactor object')}
      node <- PF$nodes[factor]
    } else {stop('PhyloFactor object contains no nodes or factors')}
  }

  nd <- which(PF$nodes==node)
    if (nd>1){
     atms <- atoms(PF$basis[,1:(nd-1),drop=FALSE])
      splt <- which(PF$basis[,nd]>0)
      grp <- atms[[which(unlist(lapply(atms,function(x,y){all(y %in% x)},y=splt)))]]

      grp1 <- intersect(grp,Descendants(tree,PF$nodes[nd],type='tips')[[1]])  ## this is our group
      grp2 <- setdiff(grp,grp1)                                               ## this is our complement
    } else {
      otus <- tree$tip.label

      grp <- 1:length(PF$basis[,nd])
      grp1 <- intersect(grp,Descendants(tree,PF$nodes[nd],type='tips')[[1]])  ## this is our group
      grp2 <- setdiff(grp,grp1)

    }

  output <- NULL


  output$group <- summary.group(PF,tree,taxonomy,factor = nd,grp1)
  output$complement <- summary.group(PF,tree,taxonomy,factor=nd,grp2)

  if (subtree==T){
    tr <- drop.tip(tree,setdiff(unlist(atms),grp))
    edgs <- rep(2,Nedge(tr))
    cols <- rep('black',Nedge(tr))

    edgG <- extractEdges(tr,taxa=tree$tip.label[grp1],type=3)
    edgs[edgG]=8
    cols[edgG]='red'
    if (prediction==T){
    dta <- rbind(output$group$PF.prediction,output$complement$PF.prediction)
    } else { dta <- t(clr(t(rbind(output$group$otuData,output$complement$otuData))))}
    phylo.heatmapAW(tree=tr,Y=dta,tipLabels = tipLabels,edge.width=edgs,edge.color=cols,...)
    output$subtree <- tr
  }

  return(output)

}

