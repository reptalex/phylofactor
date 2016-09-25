#' Performs Phylogenetic Principal Components Analysis
#' 
#' @export
#' @param Data Positive-valued data matrix whose rownames are tip-labels in the input \code{tree}. 
#' @param tree phylo-class object whose tip-labels cover the rownames of Data.
#' @param ncores optional number of cores for built-in parallelization. Be cautious of memory - each worker is sent a copy of the dataset and tree, producing ncores+1 copies of the dataset & tree.
#' @param output.edges Logical, whether or not to output the edges in the input \code{tree} corresponding to phylogenetic components
#' @param tol tolerance for compositional matrix. Rounding error in large datasets can lead to columns of compositional data not summing to 1. 
#' @param quiet Logical, whether or not to quiet warnings.
#' @return PhyCA object containing Data, tree, basis, and edges of the phylogeny corresponding to each split.
#' @examples 
#' 
#' library(phylofactor)
#' data("FTmicrobiome")
#' 
#' Data <- FTmicrobiome$PF$Data
#' tree <- FTmicrobiome$PF$tree
#' X <- FTmicrobiome$X
#' taxonomy <- FTmicrobiome$taxonomy
#' phytools::phylo.heatmap(tree,t(compositions::clr(t(Data))))
#' 
#' 
#' phca <- PhyCA(Data,tree,ncomponents = 2)
#' phcaPAR <- PhyCA(Data,tree,ncomponents=3,ncores=2)
#' 
#' phylofactor.visualize(phcaPAR,X=X,dimension=3)
#' 
#' phytools::phylo.heatmap(tree,t(compositions::clr(t(Data))))
#' edgelabels(text=1:3,edge = unlist(phcaPAR$edges),cex=2)


PhyCA <- function(Data,tree,ncores=NULL,ncomponents=NULL,output.edges=T,tol=1e-5,quiet=T){
  
  
  ################################### Housekeeping ######################################################
  if (any(Data==0)){
    if(!quiet){
      warning('Data matrix contains zeros. Zeros will be replaced with 0.65 pseudocounts')
    }
    Data[Data==0]=0.65
  }
  
  if (!(all(abs(colSums(Data)-1)<tol))){
    if(!quiet){
      warning('input Data was not compositional - will convert to compositional')
    }
    Data <- apply(Data,MARGIN=2,FUN=function(y) y/sum(y))
  }
  
  if (!(all((rownames(Data)) %in% tree$tip.label))){stop('Data contains rows not contained in input tree')}
  if (!all(rownames(Data)==tree$tip.label)){
    warning('Re-ordering data in phca$Data to fit tree. Use output data matrix for downstream analysis')
    Data <- Data[tree$tip.label,]
  }
  #########################################################################################################
  
  n=nrow(Data)
  if (is.null(ncomponents)){
    ncomponents=nrow(Data)-1
  }
 
  output <-  NULL
  output$tree <- tree
  output$Data <- Data
  output$edges <- vector(mode='list',length=ncomponents)
  totalvar <- sum(apply(apply(Data,MARGIN=2,FUN=function(y) log(y)-mean(log(y))),MARGIN=1,var))
  treeList <- list(tree)
  binList <- list(1:nrow(Data))
  
  
  #### Getting the Groups is parallelizable and very memory-intensive.
  #### To paralellize, we will 
  if (is.null(ncores)){
    Grps <- getGroups(tree)
  } else {
    cl <- phyloFcluster(ncores)
    LogData <- log(Data)
    parallel::clusterExport(cl,'LogData',envir=environment())
    nms=rownames(Data)
    treetips <- sapply(treeList,FUN=ape::Ntip)
    grpsizes <- sapply(treeList,FUN=function(tree,lg) ape::Nnode(phy=tree,internal.only=lg),lg=F)
    nnodes <- sum(grpsizes)
    cl_node_map <- sample(1:nnodes)  ### By randomizing, we can help clusters have a more even load.
    tree_map <- cumsum(grpsizes) # if tree_map[i-1]<Nde<=tree_map[i], then node is Nde-tree_map[i-1] in tree i.
    ix_cl <- parallel::clusterSplit(cl,cl_node_map)
  }
  
  
  for (phcas in 1:ncomponents){
    
    if (phcas>1){
        treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
        binList <- updateBinList(binList,grp)
        
        if (is.null(ncores)){
          Grps <- phylofactor::getNewGroups(tree,treeList,binList)
        } else {
          treetips <- sapply(treeList,FUN=ape::Ntip)
          grpsizes <- sapply(treeList,FUN=function(tree,lg) ape::Nnode(phy=tree,internal.only=lg),lg=F)
          nnodes <- sum(grpsizes)
          cl_node_map <- sample(1:nnodes)  ### By randomizing, we can help clusters have a more even load.
          tree_map <- cumsum(grpsizes) # if tree_map[i-1]<Nde<=tree_map[i], then node is Nde-tree_map[i-1] in tree i.
          ix_cl <- parallel::clusterSplit(cl,cl_node_map)
        }
    }
    
    if (is.null(ncores)){
      Y <- lapply(Grps,FUN = function(g,d) amalg.ILR(g,LogData=log(d)), d=Data)
      vars <- sapply(Y,var)
      winner <- which(vars==max(vars))
      if (length(winner)>1){
        if (!quiet){
          warning(paste('There was a tie at step ',phcas,'. The first entry will be chosen.'))
        }
        winner <- winner[1]
      }
      
      grp <- getLabelledGrp(winner,tree,Grps)
      v=ilrvec(Grps[[winner]],n)
      output$basis <- cbind(output$basis,v)
      output$projection <- rbind(output$projection,Y[[winner]])
      output$PercentVariance <- c(output$PercentVariance,vars[winner]/totalvar)
      if (output.edges){
        output$edges[[phcas]] <- list(getFactoredEdges(v,tree))
      }
      
    } else {                                                 #nset,tree_map,treeList,treetips,choice,smallglm=F,frmla=NULL,X=NULL,.
      Winners=parallel::clusterApply(cl,x=ix_cl,fun= function(x,tree_map,treeList,treetips,choice) findWinner(x,tree_map,treeList,treetips,choice),tree_map=tree_map,treeList=treeList,treetips=treetips,choice='phyca')
      vs <- sapply(Winners,FUN=function(x) x$var)
      winner=which(vs==max(vs))
      if (length(winner)>1){
        if (!quiet){
          warning(paste('There was a tie at step ',phcas,'. The first entry will be chosen.'))
        }
        winner <- winner[1]
      }
      
      Winners <- Winners[[winner]]
      grp <- lapply(Winners$grp,FUN=function(x,nms) which(nms %in% x),nms=nms)
      lns <- lapply(grp,length)
      names(grp)[lns==1]='tip'
      v <- phylofactor::ilrvec(grp,n)
      output$basis <- cbind(output$basis,v)
      output$projection <- rbind(output$projection,Winners$Y)
      output$PercentVariance <- c(output$PercentVariance,Winners$var/totalvar)
      if (output.edges){
        output$edges[[phcas]] <- list(phylofactor::getFactoredEdges(v,tree))
      }
    }

  }
  
  if(!is.null(ncores)){
    parallel::stopCluster(cl)
    rm('cl')
    gc()
  }
  
  class(output) <- 'phyca'
  return(output)
  
}