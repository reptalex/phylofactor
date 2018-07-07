#' PhyloFactor internal function to prep, find and summarize the best ILR coordinate
#' @export
#' @param TransformedData Transformed data matrix whose rows are species and columns are samples.
#' @param X independent variables for input into glm.
#' @param frmla Formula for input into glm by lapply(Y,FUN = pglm,x=xx,frmla=frmla,choice,...). 
#' @param Grps Groups - a list whose elements are two-element lists containing the groups and their compliments for amalgamation into ILR coordinates
#' @param contrast.fcn Function taking as input \code{contrast.fcn(Grps,TransformedData)} and returning a matrix with rows corresponding to \code{Grps} and columns corresponding to columns in \code{TransformedData}. Must work for both list of \code{Grps} and a single element, e.g. \code{Grps[[1]]}. For advanced phylofactorization with customized \code{choice.fcn}, the output from \code{contrast.fcn} is input directly into \code{choice.fcn}.
#' @param choice Choice function for determining the group maximizing the objective function. Currently the only allowable inputs are 'var' - minimize residual varaince - and 'F' - minimize test-statistic from anova.
#' @param treeList List of trees formed by phylofactorization and cutting trees along factored edges.
#' @param cl \code{phyloFcluster} for built-in parallelization of phylofactorization calculations.
#' @param totalvar Total variance of dataset.
#' @param ix_cl Cluster split of nodes in treeList
#' @param treetips Number of tips in treeList for quickly identifying whether nodes correspond to a root
#' @param grpsizes Number of nodes in each tree of treeList
#' @param tree_map Cumulative number of nodes in trees of treeList - allows rapid mapping of nodes in ix_cl to appropriate tree in treeList.
#' @param nms rownames of TransformedData, allowing reliable mapping of rows of data to tree.
#' @param choice.fcn optional customized choice function to choose 'best' edge; see \code{\link{PhyloFactor}}
#' @param method See \code{\link{PhyloFactor}}
#' @param ... optional input arguments for \code{\link{glm}}
PhyloRegression <- function(TransformedData,X,frmla,Grps=NULL,contrast.fcn=NULL,choice,treeList=NULL,cl,totalvar=NULL,ix_cl,treetips=NULL,grpsizes=NULL,tree_map=NULL,nms=NULL,choice.fcn,method='glm',...){
   #cl - optional phyloCluster input for parallelization of regression across multiple groups.
  D <- dim(TransformedData)[1]
  xx=X
  output <- NULL
  ############# REGRESSION ################
  ##### SERIAL #####
  if (is.null(cl)){
    ngrps <- length(Grps)
    if (is.null(contrast.fcn)){
      V <- sapply(Grps,FUN=function(g,D) ilrvec(g,D),D=D) %>% t()
      Y <- V %*% TransformedData
    } else {
      Y <- contrast.fcn(Grps,TransformedData)
    }
    if (choice != 'custom'){
    ################ DEFAULT choice.fcn #########################
      if (method != 'max.var'){
        GLMs <- apply(Y,1,FUN = pglm,xx=X,frmla=frmla,...)
        stats <- matrix(sapply(GLMs,FUN=phylofactor::getStats),ncol=3,byrow=T) 
        #contains Pvalues, F statistics, and explained var
        rownames(stats) <- 1:ngrps
        colnames(stats) <- c('Pval','F','ExplainedVar')
        Yhat <- lapply(GLMs,stats::predict)
      } else {
        objective <- apply(Y,1,stats::var)
        stopStatistics <- objective
      }
    ###############################################################
    } else {
    ##################### input choice.fcn ########################
      ch <- apply(Y,1,FUN=function(y,xx,...) choice.fcn(y=y,X=xx,...),xx=xx)
      objective <- sapply(ch,getElement,'objective')
      stopStatistics <- sapply(ch,getElement,'stopStatistics')
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
    Winners=parallel::clusterApply(cl,x=ix_cl,fun= function(x,tree_map,treeList,treetips,contrast.fcn,choice,method,frmla,xx,choice.fcn,...) findWinner(x,tree_map=tree_map,treeList=treeList,treetips=treetips,contrast.fcn=contrast.fcn,choice=choice,method=method,frmla=frmla,xx=xx,choice.fcn=choice.fcn,...) ,tree_map=tree_map,treeList=treeList,treetips=treetips,contrast.fcn=contrast.fcn,choice=choice,method=method,frmla=frmla,xx=xx,choice.fcn=choice.fcn,...)

    grps <- lapply(Winners,getElement,'grp')
    Y <- lapply(grps,BalanceContrast,TransformedData=TransformedData)
    
    
    ####################################### DEFAULT REGRESSIONS #####################
    if (choice != 'custom'){
      if (method != 'max.var'){
        gg <- lapply(Y,FUN = pglm,xx=X,frmla=frmla,...)
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
      objective <- parallel::parLapply(cl,Y,choice.fcn,X,...) %>% sapply(getElement,'objective')
    #################################################################################
    }
    
    winner=which(objective==max(objective))
    if (length(winner)>1){
      warning(paste('Objective function produced',toString(length(winner)),
                      'identical groups. Will choose group at random.',sep=' '))
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
        output$stopStatistics <- stats[,'Pval']   #this can allow us to do a KS test on P-values as a stopping function for PhyloFactor
        output$explainedvar <- stats[winner,'ExplainedVar']/totalvar
      } else {
        output$explainedvar <- objective[winner]/totalvar
        output$stopStatistics <- stopStatistics
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
        output$stopStatistics <- unlist(c(sapply(Winners,FUN=function(x) x$stopStatistics)))
      }
      if (choice=='var'){
        if (!method=='max.var'){
          output$explainedvar <- objective[winner]/totalvar
        } else {
          output$explainedvar <- stats::var(Y[[winner]])/totalvar
        }
      }
      ###############################
    
    } else {
      ######## choice.fcn input #####
      for (i in 1:length(Winners)){
        output$stopStatistics <- c(output$stopStatistics,Winners[[i]]$stopStatistics)
      }
        # output$stopStatistics <- sapply(Winners,getElement,'stopStatistics')
        output$custom.output <- choice.fcn(y=Y[[winner]],X=xx,PF.output=T,...)
        ###############################
    }
    
    output$grp <- lapply(grps[[winner]],FUN=function(x,nms) which(nms %in% x),nms=nms)
    output$basis <- ilrvec(output$grp,D)
    
  }


  return(output)
}
