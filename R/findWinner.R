#' Internal phylofactor function for finding the winning edge.
#'
#' @export
#' @param nset set of nodes
#' @param tree_map mapping cumulative number of nodes in treeList, used to map elements of nset to their appropriate tree in treeList.
#' @param treeList list containing disjoint trees from phylofactor / PhyCA
#' @param treetips number of tips in each tree
#' @param LogData logarithm of data - taking logarithm beforehand allows us to compute the logarithm of big datasets only once. 
#' @param choice string indicating how we choose the winner. Must be either \code{'var'}, \code{'F'}, or \code{'phyca'}
#' @param smallglm Logical - whether or not to use regular GLM. if smallglm=F, will use bigglm from the biglm package.
findWinner <- function(nset,tree_map,treeList,treetips,choice,smallglm=F,frmla=NULL,xx=NULL,...){
  
  
  ########### set-up and prime variables #############
  grp <- vector(mode='list',length=2)
  output <- NULL
  Y <- numeric(ncol(LogData))
  if (!exists('gg')){
    gg <- NULL #This will be our GLM
  }
  
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
      if ((nn-tree_map[[whichTree-1]])==treetips[whichTree]+1){ 
        #This prevents us from drawing the root of a subtree, which has no meaningful ILR transform. 
        next 
      }
      grp[[1]] <- phangorn::Descendants(treeList[[whichTree]],node=(nn-tree_map[whichTree-1]))[[1]]
    } else {
      whichTree <- 1
      if (nn==(treetips[1]+1)){ 
        #This prevents us from drawing the root of a subtree, which has no meaningful ILR transform. 
        next 
      }
      grp[[1]] <- phangorn::Descendants(treeList[[1]],node=nn)[[1]]
    }
    
    
      grp[[2]] <- setdiff(1:treetips[whichTree],grp[[1]])
      grp <- lapply(grp,FUN=function(x,tree) tree$tip.label[x],tree=treeList[[whichTree]])
      #This converts numbered grps of tip-labels for trees in treeList to otus that correspond to rownames in Data.
      
      
      #### ILR-transform the data, splitting grp ###
      r = length(grp[[1]])
      s = length(grp[[2]])
      if (r>1){
        Y <- colSums(LogData[grp[[1]],])*(sqrt(s/(r*(r+s))))
      } else {
        Y <- LogData[grp[[1]],]*(sqrt(s/(r*(r+s))))
      }
      if (s>1){
        Y <- Y-colSums(LogData[grp[[2]],])*sqrt(r/(s*(r+s)))
      } else {
        Y <- Y-LogData[grp[[2]],]*sqrt(r/(s*(r+s)))
      }
      ##################################################
      
      if (choice %in% c('var','F')){
          
          if (!exists('dataset')){
            dataset <- c(list(Y),as.list(xx))
            names(dataset) <- c('Data',names(xx))
            dataset <- model.frame(frmla,data = dataset)
          } else {
            if (nrow(dataset) != length(Y)){
              stop(paste('dataset=',toString(dim(dataset)),' Y=',toString(length(Y)),'.',sep=''))
            }
            dataset$Data <- Y
          }
        
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
            }
          } else {
            if (stats['F']>output$Fstat){
              output$grp <- grp
              output$Fstat <- stats['F']
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
  
  if (choice %in% c('var','F') && !smallglm){ #convert bigglm to glm
    Y <- amalg.ILR(output$grp,LogData=LogData)
    dataset <- c(list(Y),as.list(xx))
    names(dataset) <- c('Data',names(xx))
    dataset <- model.frame(frmla,data = dataset)
    output$glm <- glm(frmla,data = dataset,...)
  }
  return(output)
}