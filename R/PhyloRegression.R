#' Performs regression of groups of Data against X for an input list of groups, and outputs the best group based on "choice" function input.
#' @export
#' @param Data Compositional data matrix whose rows are parts and columns are samples.
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

PhyloRegression <- function(Data,X,frmla,Grps,choice,cl,totalvar=NULL,...){
   #cl - optional phyloCluster input for parallelization of regression across multiple groups.
  n <- dim(Data)[1]
  ngrps <- length(Grps)
  #pre-allocate
  Y <- vector(mode='list',length=ngrps)
  GLMs <- Y
  stats <- matrix(NA,ncol=2,nrow=ngrps)
  
  ############# REGRESSION ################
  ##### SERIAL #####
  if (is.null(cl)){
    Y <- lapply(X=Grps,FUN=amalg.ILR,Log.Data=log(Data))
    # GLMs <- lapply(X=Y,FUN = pglm,xx=X,frmla=frmla,smallglm=T)
    GLMs <- lapply(X=Y,FUN = pglm,xx=X,frmla=frmla,smallglm=T,...)
    stats <- matrix(unlist(lapply(GLMs,FUN=getStats)),ncol=3,byrow=T) #contains Pvalues and F statistics
    rownames(stats) <- names(GLMs)
    colnames(stats) <- c('Pval','F','ExplainedVar')
    Yhat <- lapply(GLMs,predict)
  } else {  ##### PARALLEL #####
    
      if (length(Grps)>=(2*length(cl))){
        # dum <- phyloregPar(Grps,Data,X,frmla,cl,choice)
        dum <- phyloregPar(Grps,Data,X,frmla,cl,choice,...)
        Y <- dum$Y
        stats <- dum$stats #contains Pvalues and F statistics
        Yhat <- dum$Yhat #Contains predicted ilr coordintes, unless choice='F' - in that case, we only calculate prediction for the winner
      } else { #If we don't have many groups, there's no major need to parallelize. To avoid the hassle, I just serialize the computation for few groups.
        Y <- lapply(X=Grps,FUN=amalg.ILR,Log.Data=log(Data))
        GLMs <- lapply(X=Y,FUN = pglm,xx=X,frmla=frmla,smallglm=T,...)
        stats <- matrix(unlist(lapply(GLMs,FUN=getStats)),ncol=3,byrow=T) #contains Pvalues and F statistics
        rownames(stats) <- names(GLMs)
        colnames(stats) <- c('Pval','F','ExplainedVar')
        Yhat <- lapply(GLMs,predict)
      }
    
  }
  

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
  



  ############ OUTPUT ##########################
  output <- NULL
  output$winner <- winner
  output$basis <- ilrvec(Grps[[winner]],n) #this allows us to quickly project other data onto our partition

  if (is.null(cl)){
    output$glm <- GLMs[[winner]]         #this will enable us to easily extract effects and contrasts between clades, as well as project beyond our dataset for quantitative independent variables.
  } else {
    output$glm <- pglm(Y[[winner]],X,frmla,smallglm=T,...)
  }

  output$p.values <- stats[,'Pval']   #this can allow us to do a KS test on P-values as a stopping function for PhyloFactor

 if (choice=='var'){
      output$explainedvar <- stats[winner,'ExplainedVar']/totalvar
 }

  output$residualData <- PredictAmalgam(predict(pglm(Y[[winner]],xx=X,frmla=frmla,smallglm=T)),Grps[[winner]],n)

  return(output)
}
