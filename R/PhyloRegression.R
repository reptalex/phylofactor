#' Performs regression of groups of Data against X for an input list of groups, and outputs the best group based on "choice" function input.
#' @export
#' @param Data Compositional data matrix whose rows are parts and columns are samples.
#' @param X independent variables for input into glm.
#' @param frmla Formula for input into glm by lapply(Y,FUN = pglm,x=X,frmla=frmla,choice,...)
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

PhyloRegression <- function(Data,X,frmla,Grps,choice,cl,Pval.Cutoff,...){
  #Groups taxa by "Grps" and regresses independent variable X on Data according to formula.
  #Regression is done for a set of groups, "Grps", out of a total possible set of taxa, "set" (e.g. can perform regression on log-ratio (1,2) over (3,4) out of (1,2,3,4,5))
  #Data - data matrix - rows are otus and columns correspond to X.
  #X - independent variable
  #frmla - object of class "formula" indicating the regression of variable "Data" in terms of variable "X".
  #Grps - set of groups (must be list of 2-element lists containing taxa within 'set')
     #'add' looks at the log-ratio of relative abundances of the taxa in the two groups of Grps[]
     #'ILR' uses geometric means as centers of groups, and regression is performed on balances of groups according to the ILR method of Egozcue et al. (2003)
  #choice - method for choosing the dominant partition in tree, either 't' or 'var'.
     #'t' will choose dominant partition based on the Grps whose regression maximized the test-statistic
     #'var' will choose based on Grps which maximized the percent explained variance in the clr-transformed dataset.
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
    stats <- matrix(unlist(lapply(GLMs,FUN=getStats)),ncol=2,byrow=T) #contains Pvalues and F statistics
    rownames(stats) <- names(GLMs)
    colnames(stats) <- c('Pval','F')
    Yhat <- lapply(GLMs,predict)
  } else {  ##### PARALLEL #####
    
      if (length(Grps)>=(2*length(cl))){
        # dum <- phyloregPar(Grps,Data,X,frmla,cl,Pval.Cutoff)
        dum <- phyloregPar(Grps,Data,X,frmla,cl,choice,Pval.Cutoff,...)
        # GLMs <- dum$GLMs
        Y <- dum$Y
        stats <- dum$stats #contains Pvalues and F statistics
        Yhat <- dum$Yhat #Contains predicted ilr coordintes, unless choice='F' - in that case, we only calculate prediction for the winner
      } else { #If we don't have many groups, there's no major need to parallelize. To avoid the hassle, I just serialize the computation for few groups.
        Y <- lapply(X=Grps,FUN=amalg.ILR,Log.Data=log(Data))
        GLMs <- lapply(X=Y,FUN = pglm,xx=X,frmla=frmla,smallglm=T,...)
        stats <- matrix(unlist(lapply(GLMs,FUN=getStats)),ncol=2,byrow=T) #contains Pvalues and F statistics
        rownames(stats) <- names(GLMs)
        colnames(stats) <- c('Pval','F')
        Yhat <- lapply(GLMs,predict)
      }
  }
  

  ############# CHOICE - EXTRACT "BEST" CLADE ##################
  any.sigs=T  #This is used to track the P-value cutoff. If there are no Pvalues < Pval.Cutoff, any.sigs=F below.
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
      
      #calculate only the residual vars for those Pvals <= Pval.Cutoff
      if (is.null(Pval.Cutoff)==F){
        Ps <- which(stats[,'Pval']<=Pval.Cutoff)
        if (length(Ps)==0){
          any.sigs = F
        }
      } else {
        Ps <- 1:ngrps
      }
      
      if (is.null(cl)){
        
        #Pre-allocate - in the future, we're going to remember the residual variance for groups from previous, un-split trees/bins.
        if (!exists('residualvar')){
          residualvar <- rep(Inf,ngrps)
        }
        if (!exists('predictions') && length(Ps)>0){
          predictions <- vector(mode='list',length=length(Ps))

        }
        
        
        #### Calculate predicted dataset and residual variance 
        if (length(Ps)>0){
        predictions <- mapply(PredictAmalgam,Yhat[Ps],Grps[Ps],n,SIMPLIFY=F)
        residualvar[Ps] <- sapply(predictions,residualVar,Data=Data)
        }
        winner <- which(residualvar == min(residualvar))
      } else {
        if (length(Grps)>=(2*length(cl))){
          winner <- which(dum$residualvar==min(dum$residualvar))
        } else {
          
          if (!exists(residualvar)){
            residualvar <- rep(Inf,ngrps)
          }
          if (!exists(predictions)){
            predictions <- vector(mode='list',length=length(Ps))
          }
          
          predictions <- mapply(PredictAmalgam,Yhat[Ps],Grps[Ps],n,SIMPLIFY=F)
          residualvar[Ps] <- sapply(predictions,residualVar,Data=Data)
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
  
  if (!any.sigs){
    output$STOP=T
  } else {
    output$STOP=F
  }
  
 
    output$basis <- ilrvec(Grps[[winner]],n) #this allows us to quickly project other data onto our partition


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

  output$residualData <- PredictAmalgam(predict(pglm(Y[[winner]],xx=X,frmla=frmla,smallglm=T)),Grps[[winner]],n)

  return(output)
}
