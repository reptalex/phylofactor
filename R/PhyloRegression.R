#' Performs regression of groups of Data against X for an input list of groups, and outputs the best group based on "choice" function input.
#' @export
#' @param Data Compositional data matrix whose rows are parts and columns are samples.
#' @param X independent variables for input into glm.
#' @param frmla Formula for input into glm by lapply(Y,FUN = pglm,x=X,frmla=frmla,choice,...)
#' @param Grps Groups - a list whose elements are two-element lists containing the groups and their compliments for log-ratio regression by "method"
#' @param method Method for amalgamation and comparison of groups. Default is method='ILR'.
#' @param choice Choice function for determining the group maximizing the objective function. Currently the only allowable inputs are 'var' - minimize residual varaince - and 'F' - minimize test-statistic from anova.
#' @param cl phyloFcluster input for built-in parallelization of grouping, amalgamation, regression, and objective-function calculation.
#' @param Pbasis Coming soon - input Pbasis for amalgamation method "add".
#' @examples
#' data("FTmicrobiome")
#' Y <- FTmicrobiome$OTUTable
#' Y <- Y[which(rowSums(Y==0)<30),] #only include taxa present in at least 29 samples
#' Y[Y==0]=.65
#'
#' Y <- t(clo(t(Y)))
#' tree <- drop.tip(FTmicrobiome$tree, setdiff(tree$tip.label,rownames(Y)))
#' X <- FTmicrobiome$X
#' Grps <- getGroups(tree)
#' method='ILR'
#' choice='var'
#' cl <- phyloFcluster(2)
#' pr <- PhyloRegression(Y,X,frmla,Grps,method,choice,cl)
#'
#' stopCluster(cl)
#' gc()
#' par(mfrow=c(1,2))
#' image(clr(t(Y)),main="Original Data")
#' image(clr(t(pr$residualData)),main="residual Data")

PhyloRegression <- function(Data,X,frmla,Grps,method,choice,cl,Pbasis=1,...){
  #Groups taxa by "Grps" and regresses independent variable X on Data according to formula.
  #Regression is done for a set of groups, "Grps", out of a total possible set of taxa, "set" (e.g. can perform regression on log-ratio (1,2) over (3,4) out of (1,2,3,4,5))
  #Data - data matrix - rows are otus and columns correspond to X.
  #X - independent variable
  #frmla - object of class "formula" indicating the regression of variable "Data" in terms of variable "X".
  #Grps - set of groups (must be list of 2-element lists containing taxa within 'set')
  #method - method for amalgamation of groups, either 'add' or 'multiply'.
     #'add' looks at the log-ratio of relative abundances of the taxa in the two groups of Grps[]
     #'ILR' uses geometric means as centers of groups, and regression is performed on balances of groups according to the ILR method of Egozcue et al. (2003)
  #choice - method for choosing the dominant partition in tree, either 't' or 'var'.
     #'t' will choose dominant partition based on the Grps whose regression maximized the test-statistic
     #'var' will choose based on Grps which maximized the percent explained variance in the clr-transformed dataset.
  #cl - optional phyloCluster input for parallelization of regression across multiple groups.
  if(is.null(Pbasis)){Pbasis=1}
  n <- dim(Data)[1]

  ############# REGRESSION ################
  #these can both be parallelized
  if (is.null(cl)){
    Y <- lapply(X=Grps,FUN=amalgamate,Data=Data,method)
    GLMs <- lapply(X=Y,FUN = pglm,x=X,frmla=frmla,smallglm=T,...)
    stats <- matrix(unlist(lapply(GLMs,FUN=getStats)),ncol=2,byrow=T) #contains Pvalues and F statistics
    rownames(stats) <- names(GLMs)
    colnames(stats) <- c('Pval','F')
    Yhat <- lapply(GLMs,predict)
  } else {

    if (length(Grps)>=(2*length(cl))){
      ## the following includes paralellization of residual variance if choice=='var'
      # dum <- phyloregPar(Grps,Data,X,frmla,choice,method,Pbasis,cl,...)
      dum <- phyloregPar(Grps,Data,X,frmla,choice,method,Pbasis,cl)
      # GLMs <- dum$GLMs
      Y <- dum$Y
      stats <- dum$stats #contains Pvalues and F statistics
      Yhat <- dum$Yhat
    } else {
      Y <- lapply(X=Grps,FUN=amalgamate,Data=Data,method)
      GLMs <- lapply(X=Y,FUN = pglm,x=X,frmla=frmla,smallglm=T,...)
      stats <- matrix(unlist(lapply(GLMs,FUN=getStats)),ncol=2,byrow=T) #contains Pvalues and F statistics
      rownames(stats) <- names(GLMs)
      colnames(stats) <- c('Pval','F')
      Yhat <- lapply(GLMs,predict)
    }

  }

  ############# CHOICE - EXTRACT "BEST" CLADE ##################
  if (choice=='F'){
    #in this case, the GLM outputs was an F-statistic.
    winner <- which(stats[,'F']==max(stats[,'F']))
  } else { #we pick the clade that best reduces the residual variance. Since all have the same total variance
           #this means we pick the clade with the minimum residual variance.
           #To be clear, here for total variance of a data matrix, X, I'm using var(c(X)).

      #the  GLM.outputs above are predictions on log-ratios based on the "amalgamation" function.
      #The computation of percent variance explained by each
      #partition can also be parallelized. A function such as
      #phyloregPar and phyloregParV will be useful to parellelize
      # Y, GLM.outputs, and PercVar
      totalvar <- var(c(compositions::clr(t(Data))))
      if (is.null(cl)){
        predictions <- mapply(PredictAmalgam,Yhat,Grps,n,method,Pbasis,SIMPLIFY=F)
        residualvar <- sapply(predictions,residualVar,Data=Data)
        winner <- which(residualvar == min(residualvar))
      } else {
        if (length(Grps)>=(2*length(cl))){
          winner <- which(dum$residualvar==min(dum$residualvar))
        } else {
          predictions <- mapply(PredictAmalgam,Yhat,Grps,n,method,Pbasis,SIMPLIFY=F)
          residualvar <- sapply(predictions,residualVar,Data=Data)
          winner <- which(residualvar == min(residualvar))
        }
      }
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

  if (method=='ILR'){
    output$basis <- ilrvec(Grps[[winner]],n) #this allows us to quickly project other data onto our partition
  } else { ### need to build basis for method='add'
    output$basis <- .............
  }

  if (is.null(cl)){
    output$glm <- GLMs[[winner]]         #this will enable us to easily extract effects and contrasts between clades, as well as project beyond our dataset for quantitative independent variables.
  } else {
    output$glm <- pglm(Y[[winner]],X,frmla,smallglm=T,...)
  }

  output$p.values <- stats[,'Pval']   #this can allow us to do a KS test on P-values as a stopping function for PhyloFactor

 if (choice=='var'){
    if (is.null(cl)){
      output$explainedvar <- residualvar[winner]/totalvar
    } else {
      if (length(Grps)>=(2*length(cl))){
        output$explainedvar <- dum$residualvar[winner]/totalvar
      } else {
        output$explainedvar <- residualvar[winner]/totalvar
      }
    }
 }

  output$residualData <- PredictAmalgam(Yhat[[winner]],Grps[[winner]],n,method,Pbasis)

  return(output)
}
