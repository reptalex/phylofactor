#' Performs regression of groups of Data against X for an input list of groups, and outputs the best group based on "choice" function input.
#' @export
#' @param LogData Logarithm of Compositional data matrix whose rows are parts and columns are samples.
#' @param X independent variables for input into glm.
#' @param frmla Formula for input into glm by lapply(Y,FUN = pglm,x=X,frmla=frmla,choice,...). 
#' @param Grps Groups - a list whose elements are two-element lists containing the groups and their compliments for amalgamation into ILR coordinates
#' @param choice Choice function for determining the group maximizing the objective function. Currently the only allowable inputs are 'var' - minimize residual varaince - and 'F' - minimize test-statistic from anova.
#' @param cl phyloFcluster input for built-in parallelization of grouping, amalgamation, regression, and objective-function calculation.
#' @examples
#' data("FTmicrobiome")
#' library(ape)
#' library(compositions)
#' Y <- FTmicrobiome$OTUTable
#' Y <- Y[which(rowSums(Y==0)<30),] #only include taxa present in at least 29 samples
#' Y[Y==0]=.65
#'
#' Y <- t(clo(t(Y)))
#' tree <- drop.tip(FTmicrobiome$tree, setdiff(FTmicrobiome$tree$tip.label,rownames(Y)))
#' X <- FTmicrobiome$X
#' Grps <- getGroups(tree)
#' choice='var'
#' cl <- phyloFcluster(2)
#'
#' frmla <- Y~X
#' print(head(Y))
#' Z=Y
#'
#' pr <- PhyloRegression(Z,X,frmla,Grps,choice,cl)
#'
#' stopCluster(cl)
#' gc()
#' par(mfrow=c(1,2))
#' image(clr(t(Y)),main="Original Data")
#' image(clr(t(pr$residualData)),main="residual Data")

PhyloRegression <- function(LogData,X,frmla,Grps,choice,treeList,cl,totalvar=NULL,ix_cl,treetips,grpsizes,tree_map,quiet=T,nms=NULL,...){
   #cl - optional phyloCluster input for parallelization of regression across multiple groups.
  n <- dim(LogData)[1]
  
  
  ############# REGRESSION ################
  ##### SERIAL #####
  if (is.null(cl)){
    ngrps <- length(Grps)
    #pre-allocate
    Y <- vector(mode='list',length=ngrps)
    GLMs <- Y
    stats <- matrix(NA,ncol=2,nrow=ngrps)
    Y <- lapply(X=Grps,FUN=amalg.ILR,LogData=LogData)
    GLMs <- lapply(X=Y,FUN = pglm,xx=X,frmla=frmla,smallglm=T,...)
    stats <- matrix(unlist(lapply(GLMs,FUN=getStats)),ncol=3,byrow=T) #contains Pvalues and F statistics
    rownames(stats) <- names(GLMs)
    colnames(stats) <- c('Pval','F','ExplainedVar')
    Yhat <- lapply(GLMs,predict)
    
    
      ############# CHOICE - EXTRACT "BEST" CLADE ##################
      if (choice=='F'){ #in this case, the GLM outputs was an F-statistic.
          winner <- which(stats[,'F']==max(stats[,'F']))
      } else { #we pick the clade that best reduces the residual variance.
          winner=which(stats[,'ExplainedVar']==max(stats[,'ExplainedVar']))
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
          warning(paste(paste('minimizing residual variance produced',toString(length(winner)),sep=' '),' groups, some of which were not identical. Will choose only the first group',sep=''))
        }
        winner <- winner[1]
      }
    
    
  } else {  ##### PARALLEL #####
    # Winners=parallel::clusterApply(cl,x=ix_cl,fun= function(x,tree_map,treeList,treetips,choice,smallglm,frmla,...) findWinner(nset=x,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,smallglm=smallglm,frmla=frmla,...) ,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,smallglm=F,frmla=frmla,...)
    Winners=parallel::clusterApply(cl,x=ix_cl,fun= function(x,tree_map,treeList,treetips,choice,smallglm,frmla,xx) findWinner(nset=x,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,smallglm=smallglm,frmla=frmla,xx=xx) ,tree_map=tree_map,treeList=treeList,treetips=treetips,choice=choice,smallglm=F,frmla=frmla,xx=X)
    
    grps <- lapply(Winners,FUN=function(x) x$grp)
    Y <- lapply(grps,amalg.ILR,LogData=LogData)
    getGLMs <- function(y,frmla,X,...){
      dataset <- c(list(y),as.list(X))
      names(dataset) <- c('Data',names(X))
      dataset <- model.frame(frmla,data = dataset)
      gg=glm(frmla,data = dataset,...)
    }
    gg <- lapply(Y,FUN=getGLMs,frmla=frmla,X,...)
    stats <- lapply(gg,getStats)
    
    if (choice=='var'){
      objective <- sapply(stats,function(x) x['ExplainedVar'])
    }
    if (choice=='F'){
      objective <- sapply(stats,function(x) x['F'])
    }
    winner=which(objective==max(objective))
    
    if (length(winner)>1){
      if (!quiet){
        warning(paste('There was a tie at step ',phcas,'. The first entry will be chosen.'))
      }
      winner <- winner[1]
    }
  }

  


  ############ OUTPUT ##########################
  output <- NULL

  if (is.null(cl)){
    output$glm <- GLMs[[winner]]         #this will enable us to easily extract effects and contrasts between clades, as well as project beyond our dataset for quantitative independent variables.
    output$winner <- winner
    output$basis <- ilrvec(Grps[[winner]],n) #this allows us to quickly project other data onto our partition
    output$p.values <- stats[,'Pval']   #this can allow us to do a KS test on P-values as a stopping function for PhyloFactor
    if (choice=='var'){
      output$explainedvar <- stats[winner,'ExplainedVar']/totalvar
    }
  } else {
    output$glm <- gg[[winner]]
    output$winner <- NA
    grp <- lapply(grps[[winner]],FUN=function(x,nms) which(nms %in% x),nms=nms)
    output$grp <- grp
    output$basis <- ilrvec(grp,n)
    output$p.values <- c(sapply(Winners,FUN=function(x) x$p.values))
    if (choice=='var'){
      output$explainedvar <- objective[winner]/totalvar
    }
  }


  return(output)
}
