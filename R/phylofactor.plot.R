#' Label and visualize phylofactors
#'
#'  @export
#'  @param PF PhyloFactor object
#'  @param tree Phylogeny used in generating PhyloFactor object. Default will scan PF for a tree.
#'  @param Data Data used to generate PhyloFactor object. If Null, default will plot the phylofactor prediction of the dataset.
#'  @param bg Background color for node labels, default 'white'
#'  @param cex Size of node labels, default 2.
#'  @param clades Clades, also referred to as the splits or factors, of phylofactor object to be labelled on tree and used in prediction.
#'  @param compare Whether or not to compare the Data with phylofactor prediction. If compare = TRUE, two phylo-heatmaps will be produced for direct comparison.
#'  @param ... additional arguments for phylo.heatmap
#'  @return phylo.heatmap labelling the tree with phylofactors and optionally comparing data to phylofactor predictions.
#'  @examples
#'  ### Create Data ###
#'  set.seed(1)
#' tree <- unroot(rtree(7))

#' X <- as.factor(c(rep(0,5),rep(1,5)))
#' sigClades <- Descendants(tree,c(9,12),type='tips')
#' Data <- matrix(rlnorm(70,meanlog = 8,sdlog = .5),nrow=7)
#' rownames(Data) <- tree$tip.label
#' colnames(Data) <- X
#' Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
#' Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
#' Data <- t(clo(t(Data)))

#' PF <- PhyloFactor(Data,tree,X,nclades=2)

#' plot.phylofactor(PF,clades=1)
#' plot.phylofactor(PF,clades=c(1,2),compare=T)


phylofactor.plot <- function(PF,tree=NULL,Data=NULL,bg='white',cex=2,clades=1,compare=F,...){
  #returns phylo.heatmap highlighting our PFs.
  if(is.null(tree)){
    if(is.null(PF$tree)){stop('Input Phylofactor object does not contain tree - must input tree')}
    tree <- PF$tree
  }

  if (compare==F){
    #makes just one plot
    par(mfrow=c(1,1))
    if (is.null(Data)==T){
      if (is.null(names)){stop('must input rownames for PhyloPF predicted dataset')}
      row.names=tree$tip.label
      PData <- phylofactor.predict(PF,clades)
      rownames(PData) <- row.names
      colnames(PData) <- colnames(PData)
      phytools::phylo.heatmap(tree,t(compositions::clr(t(PData))),...)
    } else {
      phytools::phylo.heatmap(tree,t(compositions::clr(t(Data))),...)
    }
    # ape::nodelabels(text = as.list(clades),node=PF$nodes[clades],bg = bg,cex=cex)
  } else {


    #makes two plots for comparison
    par(mfrow=c(2,1))
    if (is.null(Data)==T){
      stop('if compare==T, need to input Data for Comparison')
    }
    row.names=rownames(Data)
    PData <-  phylofactor.predict(PF,clades)
    rownames(PData) <- row.names
    rownames(PData) <- tree$tip.label
    colnames(PData) <- colnames(Data)

    phytools::phylo.heatmap(tree,t(compositions::clr(t(Data))))
    # ape::nodelabels(text = as.list(clades),node=PF$nodes[clades],bg = bg,cex=cex)
    phytools::phylo.heatmap(tree,t(compositions::clr(t(PData))))
    # ape::nodelabels(text = as.list(clades),node=PF$nodes[clades],bg = bg,cex=cex)




  }
}
