#' PhyloFactor internal function to prep, find and summarize the best ILR coordinate
#' @export
#' @param LogData Logarithm of Compositional data matrix whose rows are parts and columns are samples.
#' @param X independent variables for input into glm.
#' @param frmla Formula for input into glm by lapply(Y,FUN = pglm,x=xx,frmla=frmla,choice,...). 
#' @param Grps Groups - a list whose elements are two-element lists containing the groups and their compliments for amalgamation into ILR coordinates
#' @param choice Choice function for determining the group maximizing the objective function. Currently the only allowable inputs are 'var' - minimize residual varaince - and 'F' - minimize test-statistic from anova.
#' @param treeList List of trees formed by phylofactorization and cutting trees along factored edges.
#' @param cl phyloFcluster input for built-in parallelization of grouping, amalgamation, regression, and objective-function calculation.
#' @param totalvar Total variance of dataset.
#' @param ix_cl Cluster split of nodes in treeList
#' @param treetips Number of tips in treeList for quickly identifying whether nodes correspond to a root
#' @param grpsizes Number of nodes in each tree of treeList
#' @param tree_map Cumulative number of nodes in trees of treeList - allows rapid mapping of nodes in ix_cl to appropriate tree in treeList.
#' @param quiet Logical to supress warnings
#' @param nms rownames of LogData, allowing reliable mapping of rows of data to tree.
#' @param smallglm Logical. See \code{\link{PhyloFactor}}
#' @param choice.fcn optional customized choice function to choose 'best' edge; see \code{\link{PhyloFactor}}
#' @param method See \code{\link{PhyloFactor}}
#' @param ... optional input arguments for \code{\link{glm}}
PhyloRegression <- function(LogData,X,frmla,Grps=NULL,choice,treeList=NULL,cl,totalvar=NULL,ix_cl,treetips=NULL,grpsizes=NULL,tree_map=NULL,quiet=T,nms=NULL,smallglm=F,choice.fcn,method='glm',...){
   #cl - optional phyloCluster input for parallelization of regression across multiple groups.
  D <- dim(LogData)[1]
  xx=X
  
  ############# REGRESSION ################
  ##### SERIAL #####
  if (is.null(cl)){
    ngrps <- length(Grps)
    V <- sapply(Grps,FUN=function(g,D) ilrvec(g,D),D=D) %>% t()
    Y <- V %*% LogData
    if (choice != 'custom'){
    ################ DEFAULT choice.fcn #########################
      # GLMs <- apply(Y,1,FUN = phylofactor::pglm,xx=X,frmla=frmla,smallglm=T)
      if (method != 'max.var'){
        GLMs <- apply(Y,1,FUN = pglm,xx=X,frmla=frmla,smallglm=T,...)
        stats <- matrix(sapply(GLMs,FUN=phylofactor::getStats),ncol=3,byrow=T) 
        #contains Pvalues, F statistics, and explained var
        rownames(stats) <- 1:ngrps
        colnames(stats) <- c('Pval','F','ExplainedVar')
        Yhat <- lapply(GLMs,stats::predict)
      } else {
        objective <- apply(Y,1,stats::var)
      }
    ###############################################################
    } else {
    ##################### input choice.fcn ########################
      ch <- apply(Y,1,FUN=function(y,xx,...) choice.fcn(y=y,X=xx,...),xx=xx)
      objective <- sapply(ch,FUN=function(x) x$objective)
    ###############################################################
    }
    
    ############# CHOICE - EXTRACT "BEST" CLADE ##################
    if (choice=='F'){ #in this case, the GLM outputs was an F-statistic.
      winner <- which(stats[,'F']==max(stats[,'F']))
    } else if (choice=='var' & method!='max.var'){ #we pick the clade that best reduces the residual variance.
      winner=which(stats[,'ExplainedVar']==max(stats[,'ExplainedVar']))
    } else {
      winner=which(objective==max(objective))
    }
    if (length(winner)>1){
      ## Check to see if they're equivalent groups, in which case don't display message.
      NoWarning=T
      for (nn in 1:(length(winner)-1)){
        for (mm in 2:length(winner)){
          NoWarning = NoWarning & all(Grps[nn][[1]] %in% Grps[mm][[1]]) & all(Grps[mm][[1]] %in% Grps[nn][[1]])
        }
      }
      if (NoWarning==F){
        warning(paste('Objective function produced',toString(length(winner)),
                      'identical groups. Will choose group at random.',sep=' '))
      }
      winner <- sample(winner,1)
    }
  } else {  ##### PARALLEL #####
    Winners=parallel::clusterApply(cl,x=ix_cl,fun= function(x,tree_map,treeList,treetips,choice,method,smallglm,frmla,xx,choice.fcn,...) findWinner(x,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,method=method,smallglm=smallglm,frmla=frmla,xx=xx,choice.fcn=choice.fcn,...) ,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,method=method,smallglm=smallglm,frmla=frmla,xx=xx,choice.fcn=choice.fcn,...)
    # Winners=parallel::clusterApply(cl,x=ix_cl,fun= function(x,tree_map,treeList,treetips,choice,smallglm,frmla,xx,choice.fcn) findWinner(x,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,smallglm=smallglm,frmla=frmla,xx=xx,choice.fcn=choice.fcn) ,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,smallglm=smallglm,frmla=frmla,xx=xx,choice.fcn=choice.fcn)
    
    # Winners=lapply(ix_cl,FUN=function(x,tree_map,treeList,treetips,choice,smallglm,frmla,xx,choice.fcn) findWinner(nset=x,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,smallglm=smallglm,frmla=frmla,xx=xx,choice.fcn=choice.fcn) ,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,smallglm=smallglm,frmla=frmla,xx=xx,choice.fcn=choice.fcn)
    # Recall: output from findWinner is $grp and then our objective function output: $objective, $Fstat, or $ExplainedVar, corresponding to choice='custom','F', and 'var', respectivley.
    
    grps <- lapply(Winners,FUN=function(x) x$grp)
    Y <- lapply(grps,amalg.ILR,LogData=LogData)
    
    
    ####################################### DEFAULT REGRESSIONS #####################
    if (choice != 'custom'){
      if (method != 'max.var'){
        gg <- lapply(Y,FUN = pglm,xx=X,frmla=frmla,smallglm=T,...)
        stats <- lapply(gg,getStats)
        if (choice=='var'){
          objective <- sapply(stats,function(x) x['ExplainedVar'])
        } else if (choice=='F'){
          objective <- sapply(stats,function(x) x['F'])
        }
      } else {
        objective <- sapply(Y,stats::var)
      }
    ###################################################################################
    } else {
    ####################################### input choice.fcn ########################
      objective <- sapply(Y,FUN=function(y,xx,...) choice.fcn(y=y,X=xx,...)$objective,xx=xx)
    #################################################################################
    }
    
    winner=which(objective==max(objective))
    if (length(winner)>1){
      if (!quiet){
        warning(paste('Objective function produced',toString(length(winner)),
                      'identical groups. Will choose group at random.',sep=' '))
      }
      winner <- sample(winner,1)
    }
  }

  
  

  ################################ OUTPUT ##########################
  output <- NULL

  if (is.null(cl)){ ############### SERIAL #########################
    
    ############ DEFAULT ##########
    if (choice != 'custom'){
      if (method != 'max.var'){
        output$model <- GLMs[[winner]] ## this will enable us to easily extract effects and contrasts between clades, as well as project beyond our dataset for quantitative independent variables.
        output$p.values <- stats[,'Pval']   #this can allow us to do a KS test on P-values as a stopping function for PhyloFactor
        output$explainedvar <- stats[winner,'ExplainedVar']/totalvar
      } else {
        output$explainedvar <- objective[winner]/totalvar
      }
      ###############################
    
    } else {
      ######## choice.fcn input #####
      output$stopStatistics <- lapply(ch,FUN=function(x) x$stopStatistics)
      output$custom.output <- choice.fcn(y=Y[winner,],X=xx,PF.output=T,...)
      ###############################
    }
    
    output$grp <- Grps[[winner]]
    output$basis <- ilrvec(output$grp,D) #this allows us to quickly project other data onto our partition
    
  } else {   #################### PARALLEL ##########################
    
    
    ######## Default ##############
    if (choice != 'custom'){
      if (!method=='max.var'){
        output$model <- gg[[winner]]
        output$p.values <- unlist(c(sapply(Winners,FUN=function(x) x$p.values)))
      }
      if (choice=='var'){
        if (!method=='max.var'){
          output$explainedvar <- objective[winner]/totalvar
        } else {
          output$explainedvar <- var(Y[[winner]])/totalvar
        }
      }
      ###############################
    
    } else {
      ######## choice.fcn input #####
      for (nn in 1:length(Winners)){
        output$stopStatistics <- c(output$stopStatistics,Winners[[nn]]$stopStatistics)
      }
        # output$stopStatistics <- lapply(Winners,FUN=function(x) x$stopStatistics)
        output$custom.output <- choice.fcn(y=Y[[winner]],X=xx,PF.output=T,...)
        ###############################
    }
    
    output$grp <- lapply(grps[[winner]],FUN=function(x,nms) which(nms %in% x),nms=nms)
    output$basis <- ilrvec(output$grp,D)
    
  }


  return(output)
}
