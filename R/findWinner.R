#' Internal PhyloRegression function for finding the winning edge.
#'
#' @export
#' @param nset set of nodes
#' @param tree_map mapping cumulative number of nodes in treeList, used to map elements of nset to their appropriate tree in treeList.
#' @param treeList list containing disjoint trees from phylofactor / PhyCA
#' @param treetips number of tips in each tree
#' @param LogData logarithm of data - taking logarithm beforehand allows us to compute the logarithm of big datasets only once. 
#' @param choice string indicating how we choose the winner. Must be either \code{'var'}, \code{'F'}, or \code{'phyca'}
#' @param smallglm Logical - whether or not to use regular \code{glm}. if smallglm=F, will use \code{\link{bigglm}} from the \code{\link{biglm}} package.
#' @param choice.fcn See \code{\link{PhyloFactor}}
findWinner <- function(nset,tree_map,treeList,treetips,choice,smallglm=F,frmla=NULL,xx=NULL,choice.fcn=NULL,...){
  
  
  ########### set-up and prime variables #############
  grp <- vector(mode='list',length=2)
  output <- NULL
  Y <- numeric(ncol(LogData))
  if (!exists('gg')){
    gg <- NULL #This will be our GLM
  }
  
  if (choice %in% c('var','F')){
    output$p.values <- numeric(length(nset))
  } else {
    output$stopStatistic <- numeric(length(nset))
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
  if (choice=='custom'){
    output$objective=-Inf
    if (is.null(choice.fcn)){
      if (is.null(choice.fcn)){
        choice.fcn <- function(y=NULL,X=NULL,PF.output=NULL){
          ch <- NULL
          ch$objective <- 1
          ch$stopStatistic <- 1
          return(ch)
        }
      }
    }
  }
  ####################################################
  
  
  
  
  #################################################### ITERATION OVER nset TO FIND WINNER #############################################################
  iteration=0
  for (nn in nset){
    iteration=iteration+1
    
    
    ############# getting grp - list of group & complement - for node in nset #####
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
    ###############################################################################
      
      
      
      ####################### ILR-transform the data ################################
      rr = length(grp[[1]])
      ss = length(grp[[2]])
      if (rr>1){
        Y <- colSums(LogData[grp[[1]],])*(sqrt(ss/(rr*(rr+ss))))
      } else {
        Y <- LogData[grp[[1]],]*(sqrt(ss/(rr*(rr+ss))))
      }
      if (ss>1){
        Y <- Y-colSums(LogData[grp[[2]],])*sqrt(rr/(ss*(rr+ss)))
      } else {
        Y <- Y-LogData[grp[[2]],]*sqrt(rr/(ss*(rr+ss)))
      }
      ################################################################################
      
      
      #################### Applying Choice Function to Y #############################
      ########### And updating output if objective > output$objective ################
      if (choice %in% c('var','F')){ ########### 2 of 3 default choice.fcns
          
          ################ Making data frame for regression #######
          if (!exists('dataset')){
            dataset <- c(list(Y),as.list(xx))
            names(dataset) <- c('Data',names(xx))
            dataset <- model.frame(frmla,data = dataset)
          } else {  #dataset already exists - we just need to update Data
            dataset$Data <- Y
          }
        #########################################################
        
        ############ Performing Regression ######################
        if(smallglm){
          gg=glm(frmla,data = dataset,...)
        } else {
          gg=biglm::bigglm(frmla,data = dataset,...)
        }
        #########################################################
        
        ############# Update output if objective is larger #######
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
        #########################################################
        
        
      } else if (choice=='phyca'){ #PhyCA
        v=var(Y)
        if (v>output$var){
          output$grp <- grp
          output$var <- v
          output$Y <- Y
        }
      } else {
        
        ################# choice.fcn ######################
        ch <- choice.fcn(y=Y,X=xx,PF.output=F,...)
        obj <- ch$objective
        output$stopStatistic[iteration] <- ch$stopStatistic
        ################# update if obj ###################
        if (obj>output$objective){
          output$grp <- grp
          output$objective <- obj
        }
      }  
    
  }
  #################################################### ITERATION OVER nset TO FIND WINNER #############################################################
  
  
  ################## modify output glm for default choices #################
  if (choice %in% c('var','F') && !smallglm){ #convert bigglm to glm
    Y <- amalg.ILR(output$grp,LogData=LogData)
    dataset <- c(list(Y),as.list(xx))
    names(dataset) <- c('Data',names(xx))
    dataset <- model.frame(frmla,data = dataset)
    output$glm <- glm(frmla,data = dataset,...)
  } else {
    output$glm <- gg
  }
  ##########################################################################
  
  return(output)
}