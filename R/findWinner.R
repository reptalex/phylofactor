#' Internal phylofactor function for finding the winning edge.
#'
#' @export
findWinner <- function(nset,tree_map,treeList,treetips,Log.Data,choice,smallglm=F,frmla=NULL,X=NULL,...){
  grp <- vector(mode='list',length=2)
  output <- NULL
  if (choice %in% c('var','F')){
    output$p.values <- numeric(length(nset))
  }
  
  if (choice=='var'){
    output$ExplainedVar=0
  }
  if (choice=='F'){
    output$Fstat=0
  }
  if (choice=='phyca'){
    output$var=0
  }
  
  
  iteration=0
  for (nn in nset){
    iteration=iteration+1
    
    if (nn>tree_map[1]){
      whichTree <- max(which(tree_map<nn))+1
    } else {
      whichTree <- 1
    }
    
    if (whichTree>1){
      if ((nn-tree_map[[whichTree-1]])==treetips[whichTree]+1){ 
      #This prevents us from drawing the root of a subtree, which has no meaningful ILR transform. 
      next 
      }
    } else {
      if (nn==treetips[whichTree]+1){ 
        #This prevents us from drawing the root of a subtree, which has no meaningful ILR transform. 
      next 
      }
    }
    
      if (whichTree>1){
        grp[[1]] <- phangorn::Descendants(treeList[[whichTree]],node=(nn-tree_map[whichTree-1]))[[1]]
      } else {
        grp[[1]] <- phangorn::Descendants(treeList[[whichTree]],node=nn)[[1]]
      }
      grp[[2]] <- setdiff(1:treetips[whichTree],grp[[1]])
      grp <- lapply(grp,FUN=function(x,tree) tree$tip.label[x],tree=treeList[[whichTree]])
      #This converts our numbered grps, which are local tip-labels for trees in treeList, to otus that correspond to rownames in Data.
      
      Y <- phylofactor::amalg.ILR(grp,Log.Data=Log.Data)
      
      if (choice %in% c('var','F')){
        
        dataset <- c(list(Y),as.list(X))
        names(dataset) <- c('Data',names(X))
        dataset <- model.frame(frmla,data = dataset)
        
        if(smallglm){
          gg=glm(frmla,data = dataset,...)
        } else {
          gg=biglm::bigglm(frmla,data = dataset,...)
        }
        
        stats=getStats(gg,y=Y)
        output$p.values[iteration] <- stats['Pval']
          if (choice=='var'){
            if (stats['ExplainedVar']>output$ExplainedVar){
              output$grp <- grp
              output$ExplainedVar <- stats['ExplainedVar']
              output$Pvalue <- stats['Pval']
              output$glm <- gg
              output$Y <- Y
              output$Yhat <- gg$fitted.values
            }
          } else {
            if (stats['F']>output$Fstat){
              output$grp <- grp
              output$Fstat <- stats['F']
              output$Pvalue <- stats['Pval']
              output$glm <- gg
              output$Y <- Y
              output$Yhat <- gg$fitted.values
            }
          }
      } else { #PhyCA
        v=var(Y)
        if (v>output$var){
          output$grp <- grp
          output$var <- v
          output$Y <- Y
        }
      }  
    
  }
  return(output)
}