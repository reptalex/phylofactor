summary.node <- function(PF,tree,taxonomy,node,subtree=F,prediction=T,tipLabels=F,...){
  #summarizes the IDs of taxa for a given node identified as important by PhyloFactor. If subtree==T, will also plot a subtree showing the taxa
  if ((node %in% PF$nodes)==F){stop('Node is not contained in PhyloFactor')}

  nd <- which(PF$nodes==node)
    if (nd>1){
      atms <- atoms(PF$basis[,1:(nd-1)])
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

