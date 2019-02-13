#' Generalized phylofactorization
#'
#' @export
#' @param Data data table containing columns of "Species", and terms in the \code{frmla}. If \code{algorithm=="mStable"}, \code{Data} must also include a column of "Sample" or, alternatively, \code{Data} can be a matrix whose rows are species and columns are samples and \code{MetaData} a data frame of meta-data with rows corresponding to columns of \code{Data} and the terms in \code{frmla} or non-phylo terms in \code{frmla.phylo}.
#' @param tree phylo class object containing all species in \code{Data}
#' @param frmla.phylo Formula used for method "phylo" and "mStable". Can use species-specific effects by multilevel factor "Species", and phylo-factors with "phylo". e.g. cbind(Successes,Failures)~x+Species*z+phylo*y has universal coefficient for x, species-specific coefficients for z, and phylo-factor contrasted coefficients for y.
#' @param frmla Formula. If \code{expfamily='binomial'}, the left hand side must have c(Successes,Failures)~. Otherwise, the variable for data is "Data", e.g. \code{Data~effort+phylo}
#' @param PartitioningVariables Character vector containing the variables in \code{frmla} to be used for phylofactorization. Objective function will be the sum of deviance from all variables listed here.
#' @param MetaData data frame or data table of meta-data containing variables found in \code{frmla}. If the \code{algorithm='mStable'}, the meta-data must contain a column "Sample" to enable aggregation of groups within samples.
#' @param nfactors integer for number of factors to find
#' @param ncores integers for number of cores to use for parallelization
#' @param model.fcn Regression function, such as glm, gam, glm.nb, gls. Must have column labelled "Deviance" in \code{\link{anova}} for default \code{objective.fcn}, \code{\link{pvDeviance}}.
#' @param objective.fcn Optional input objective function. Takes \code{model.fcn} output as its input, and returns a number which will be maximized to determine phylogenetic factor.
#' @param min.group.size Minimum group size; groups with one element having fewer members than \code{min.group.size} are automatically given objective -Inf; aggregation, \code{model.fcn} and \code{objective.fcn} are only run on groups containing as many or more members than \code{min.group.size}
#' @param algorithm Character, either "CoefContrast", "phylo", "mStable" or "mix". "CoefContrast" will partition the standardized coefficient matrix. "phylo" will produce \code{phylo} factors. "mStable" will use \code{phylo} factors and marginally-stable aggregation of groups. "mix" will use coefficient contrasts to identify the top alpha percent of edges and subsequently use the "phylo" algorithm for edge selection.
#' @param alpha Numeric between 0 and 1 (strictly greater than 0), indicating the top fraction of edges to use when \code{algorithm=='mix'}. Default is alpha=0.2 selecting top 20 percent of edges.
#' @param cluster.depends Character expression input to \code{eval(parse(text=cluster.depends))}. Evaluated in clusters to prime local environment - useful for customized \code{model.fcn} and \code{objective.fcn}
#' @param ... Additional arguments for \code{model.fcn}, e.g. for default \code{\link{glm}}, can use \code{family=binomial}, \code{weights}, etc. For \code{algorithm!='mStable'}, \code{subset} is not a valid optional argument due to \code{gpf} recursively subsetting based on phylogenetic factors. For \code{algorithm='mStable'}, \code{subset} indexes correspond to the Samples in order of \code{unique(Data$Sample)}
#' @return phylofactor object which can be manipulated with various \code{pf.} functions
#' @examples
#' library(phylofactor)
#' require(ggpubr)
#'
#' ilogit <- function(eta) 1/(1+exp(-eta))
#' set.seed(1)
#' m <- 50
#' n <- 200
#' tree <- rtree(m)
#' MetaData <- data.table('y'=rnorm(n),
#'                 'z'=rnorm(n,sd=0.5),
#'                 'Sample'=sapply(1:n,FUN=function(s) paste('Sample',s)),
#'                 key='Sample')
#' #we'll partition by 'y'.
#' binom.size=3
#' clade <- phangorn::Descendants(tree,75,'tips')[[1]]
#' clade2 <- phangorn::Descendants(tree,53,'tips')[[1]]
#'
#' ######## presence/absence dataset with affected clade #######
#' ## most species have higher P{present} with y
#' eta <- .5*MetaData$z+MetaData$y
#' p <- ilogit(eta)
#' M <- matrix(rbinom(m*n,binom.size,rep(p,times=m)),nrow=m,ncol=n,byrow=TRUE)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- MetaData$Sample
#'
#' #### the first clade decreases with y ####
#' eta1 <- .5*MetaData$z-MetaData$y
#' p1 <- ilogit(eta1)
#' for (species in clade){
#'    M[species,] <- rbinom(n,binom.size,p1)
#' }
#' #### the second clade weakly decreases with y ####
#' eta2 <- .5*MetaData$z-.3*MetaData$y
#' p2 <- ilogit(eta2)
#' for (species in clade2){
#'    M[species,] <- rbinom(n,binom.size,p2)
#' }
#'
#' #### Default algorithm: 'mix' #####
#' ### For default can input data matrix or data frame
#' ### with "Species", "Sample", and all relevant meta-data
#' DF <- matrix.to.phyloframe(M,data.name='Successes')
#' DF[,Failures:=binom.size-Successes]
#' setkey(DF,Sample)
#' DF <- DF[MetaData]
#'
#' ### DF must have "Species", "Sample", and the LHS of the formula input.
#' ### A separate data frame or data table, MetaData, can have "Sample" and the RHS of the formula.
#'
#' ### The default algorithm is "mix", which uses CoefContrast to limit the number of edges for
#' ### selection by algorithm 'phylo'. This algorithm has the greatest power but is also
#' ### computationally intensive. It's recommended that you input both frmla (used in CoefContrst)
#' ### and frmla.phylo (used in algorithm 'phylo').
#' ### For species-specific effects in algorithm 'phylo', you can use the variable "Species", e.g.
#' ### frmla.phylo=cbind(Successes,Failures)~Species*z+phylo*y. For universal/shared coefficients
#' ### for "z" across species, simply use frmla.phylo=cbind(Successes,Failures)~z+phylo*y
#' ### Since we're modelling a constant effect of z across species, we want to offset(z)
#' ### in the formula. Let's get the coefficients for that:
#'
#' model.z <- glm(cbind(Successes,Failures)~z,family=binomial,data=DF)
#' DF[,z.fit:=coef(model.z)['z']*z]
#'
#' pf <- gpf(DF,tree,frmla=cbind(Successes,Failures)~offset(z.fit)+y,
#'                   frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+phylo*y,
#'                     PartitioningVariables='y',
#'                     family=binomial,
#'                     nfactors=2,
#'                     ncores=2)
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
#'
#' ## can also use categorical variables as predictors, but notice the warning for coefficient
#' ## contrasts when PartitioningVariables are the formula's name for that variable. gpf will grep
#' ## for the PartitioningVariables in the names of coefficients - the warning below can be
#' ## ignored, but formulas whose terms grep each other may require specific PartitioningVariables.
#' DF[,y_factor:=as.factor(y>0)]
#' pf_categorical <- gpf(DF,tree,frmla=cbind(Successes,Failures)~offset(z.fit)+y_factor,
#'           frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+phylo*y_factor,
#'           PartitioningVariables='y_factor',
#'           family=binomial,
#'           nfactors=2,
#'           ncores=2)
#'
#' ### glm manipulation - can do everything except subset ###
#' ### For non-mStable, weights slow down the algorithm due to data.table indexing issues ###
#' # w <- sample(2,nrow(DF),T)
#' # pf.fancy <- gpf(DF,tree,frmla=cbind(Successes,Failures)~y,
#' #                    frmla.phylo=cbind(Successes,Failures)~phylo*y,
#' #                    PartitioningVariables='y',nfactors=2,ncores=2,
#' #                    family=binomial,weights=w,offset=DF$z.fit)
#' ###########################################################
#'
#' pf.non.offset <- gpf(DF,tree,frmla=cbind(Successes,Failures)~z+y,
#'                   frmla.phylo=cbind(Successes,Failures)~z+phylo*y,
#'                     PartitioningVariables='y',
#'                     family=binomial,
#'                     nfactors=2,
#'                     ncores=2)
#'
#'
#' ### The difference between the offset & non-offset is that the latter re-estimates the
#' ### coefficient for z in each partition, potentially introducing error & wasting degrees of
#' ### freedom in downstream factors.

#'
#' ### Algorithms "phylo", "mix", and "mStable" have a fairly high percent of the computation which
#' ### is parallelizable. The argument "ncores" performs built-in parallelization.
#'
#'
#' ### Another algorithm is "CoefContrast". For this algorithm, you need to input the frmla and
#' ## Partitioning Variables. CoefContrast is extremely fast and best done without parallelization
#' ## (as it is built off matrix multiplication).
#' pf <- gpf(DF,tree,frmla=cbind(Successes,Failures)~z+y,
#'                     PartitioningVariables='y',
#'                     algorithm='CoefContrast',
#'                     family=binomial,
#'                     nfactors=2)
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
#'
#' ####################### to partition on y, must have phylo* #########
#' ## Also, for inputting matrices into gpf for binomial glm or gam,
#' ## must input a list of matrices with "Successes" and "Failures":
#' Mat.List <- list('Successes'=M,'Failures'=binom.size-M)
#' pf <- gpf(Mat.List,tree,MetaData=MetaData,
#'           frmla.phylo=cbind(Successes,Failures)~z+phylo*y,nfactors=2,
#'           family=binomial,algorithm='mStable')
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
#'
#'
#' observed <- pf.heatmap(tree=tree,Data=M[,order(MetaData$y)])
#' predicted <- pf.heatmap(tree=tree,Data=ilogit(pf.predict(pf)[,order(MetaData$y)]))
#' ggarrange(observed,predicted,nrow=2)
#'
#' ################# Poisson Regression
#' set.seed(1)
#' eta <- 2*MetaData$z+MetaData$y
#' lambda <- exp(eta)
#' M <- matrix(rpois(m*n,rep(lambda,times=m)),nrow=m,ncol=n,byrow=TRUE)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- MetaData$Sample
#'
#' #### the first clade decreases with y ####
#' eta1 <- .3*MetaData$z-MetaData$y
#' lambda1 <- exp(eta1)
#' for (species in clade){
#'    M[species,] <- rpois(n,lambda1)
#' }
#' #### the second clade strongly increases with y ####
#' eta2 <- .3*MetaData$z-MetaData$y
#' lambda2 <- exp(eta2)
#' for (species in clade2){
#'    M[species,] <- rpois(n,lambda2)
#' }
#'
#'
#'
#' ##For non-binomial, use "Data" as response variable #########
#' pf <- gpf(M,tree,MetaData=MetaData,frmla.phylo=Data~phylo*(z+y),nfactors=2,family=poisson,
#'           PartitioningVariables='y',algorithm='mStable')
#' pf$factors
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
#'
#' par(mfrow=c(2,1))
#' observed <- pf.heatmap(tree=tree,Data=M[,order(MetaData$y)])
#' predicted <- pf.heatmap(tree=tree,Data=exp(pf.predict(pf)[,order(MetaData$y)]))
#' ggarrange(observed,predicted,nrow=2)
#' ### partition vector of data controlling for sample effort ###
#' set.seed(1)
#' effort <- rnorm(50)
#' eta <- effort-3
#' eta[clade] <- eta[clade]+6
#' eta[clade2] <- eta[clade2]+8
#' Data <- data.table('Species'=tree$tip.label,effort,
#'                    'Successes'=rbinom(50,1,ilogit(eta)),
#'                    'Sample'=1)
#' Data$Failures <- 1-Data$Successes
#' pf <- gpf(Data,tree,frmla.phylo=cbind(Successes,Failures)~effort+phylo,
#'           nfactors=2,algorithm='phylo',family=binomial)
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
#'
#'
#'
#' ############# generalized additive modelling ################
#' set.seed(1)
#' m <- 50
#' n <- 20
#' tree <- rtree(m)
#' MetaData <- data.table('y'=rnorm(n),
#'                        'z'=rnorm(n,sd=0.5),
#'                        'Sample'=sapply(1:n,FUN=function(s) paste('Sample',s)),
#'                        key='Sample')
#' #we'll partition by 'y'.
#' binom.size=10
#' clade <- phangorn::Descendants(tree,75,'tips')[[1]]
#' clade2 <- phangorn::Descendants(tree,53,'tips')[[1]]
#'
#' ######## presence/absence dataset with affected clade #######
#' ## most species have higher P{present} with y
#' eta <- .5*MetaData$z+MetaData$y
#' p <- ilogit(eta)
#' M <- matrix(rbinom(m*n,binom.size,rep(p,times=m)),nrow=m,ncol=n,byrow=TRUE)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- MetaData$Sample
#'
#' #### the first clade decreases with y ####
#' eta1 <- .5*MetaData$z+2*MetaData$y^2         #non-linear response
#' p1 <- ilogit(eta1)
#' for (species in clade){
#'   M[species,] <- rbinom(n,binom.size,p1)
#' }
#' #### the second clade weakly decreases with y ####
#' eta2 <- .5*MetaData$z-4*MetaData$y           #linear response
#' p2 <- ilogit(eta2)
#' for (species in clade2){
#'   M[species,] <- rbinom(n,binom.size,p2)
#' }
#'
#' DF <- matrix.to.phyloframe(M,data.name='Successes')
#' DF[,Failures:=binom.size-Successes]
#' setkey(DF,Sample)
#' DF <- DF[MetaData]
#'
#' model.z <- glm(cbind(Successes,Failures)~z,family=binomial,data=DF)
#' DF[,z.fit:=coef(model.z)['z']*z]
#'
#' pf <- gpf(DF,tree,frmla.phylo=cbind(Successes,Failures)~offset(z.fit)+s(y,by=phylo),
#'           PartitioningVariables='y',family=binomial,
#'           nfactors=2,ncores=2,model.fcn = mgcv::gam,algorithm = 'phylo')
#' pf$factors
gpf <- function(Data,tree,frmla.phylo=NULL,frmla=NULL,PartitioningVariables=NULL,MetaData=NULL,nfactors=NULL,ncores=NULL,model.fcn=stats::glm,objective.fcn=pvDeviance,min.group.size=1,algorithm='mix',alpha=0.2,cluster.depends='',...){

  output <- NULL
  output$call <- match.call()
  output$model.fcn <- model.fcn
  output$additional.arguments <- list(...)
  output$algorithm <- algorithm

# Check Formulas --------------------------------------------------------
  ########## checking formulas ############
  # if (is.null(frmla) & is.null(frmla.phylo)){
  #   frmla.phylo <- Data~phylo
  # }

  ### Coefficient contrast requires frmla and partitioning variables.
  ### The code below parses out input frmla.phylo into partitioning variables
  ### and frmla
  if (is.null(frmla) & algorithm %in% c('CoefContrast','mix')){
    ## We'll try to automatically parse the frmla.phylo as defining partitioning variables for CoefContrast.
    RHS <- unlist(sapply(as.character(frmla.phylo[[3]]),strsplit,'\\+'))
    RHS <- RHS[RHS!='*']
    names(RHS) <- NULL
    RHS <- RHS[!RHS=='']
    RHS <- RHS[!grepl(', by = phylo',RHS)]
    RHS <- gsub(' ','',RHS) %>% gsub('\\*',':',.)
    trms <- unique(c(base::attr(stats::terms(frmla.phylo),'term.labels'),RHS))

    if (!any(grepl('phylo',trms))){
      stop(paste('For algorithm=',algorithm,', must either input frmla and (optional) partitioning variables or frmla.phylo with phylo* or *phylo term indicating the partitioning variables',sep=''))
    } else { ## we'll define partitioning variables by phylo*, *phylo, or s(var,by=phylo)
      pure.phylo <- trms=='phylo'
      phylo.product <- grepl(':',trms) & grepl('phylo',trms)
      if (any(phylo.product)){
        product.vars <- gsub('phylo','',gsub(':','',trms[phylo.product]))
        if (length(product.vars)==0){
          stop('could not parse out phylo* or *phylo terms in frmla.phylo')
        } else {
          warning('Explanatory variables {',paste(product.vars,collapse=','),'} defined as partitioning variables for regression coefficient contrasts')
        }
        if (is.null(PartitioningVariables)){
          PartitioningVariables <- product.vars
        } else {
          if (!all(sort(PartitioningVariables)==sort(product.vars))){
            stop('Partitioning Variables and frmla.phylo do not match up for algorithm="CoefContrast" or "mix". All terms with phylo* or *phylo must be partitioning variables for coefficient contrasts.')
          }
        }
        ### must include: offset, non-partitioning & partitioning terms
        RHS <- paste(trms[!grepl('phylo',trms)],collapse='+')
        LHS <- as.character(frmla.phylo[[2]])
        if (length(LHS)==1){
          frmla <- stats::as.formula(paste(LHS,'~',RHS,sep=''))
        } else {
          if (!all.equal(c('cbind','Successes','Failures'),LHS)){
            stop('unknown left-hand side of frmla.phylo. If not a single response variable, it must be cbind(Successes,Failures) for binomial regression.')
          } else {
            frmla <- stats::as.formula(paste('cbind(Successes,Failures)~',RHS,sep=''))
          }
        }
      }

      phylo.smooth <- grepl(', by = phylo',trms)
      if (any(phylo.smooth)){
        smooth.vars <- gsub(', by = phylo','',trms[phylo.smooth]) %>%
                      strsplit(.,'\\(') %>%
                      lapply(.,FUN=function(a) gsub(')','',a)) %>%
                      sapply('[',2)
        if (length(smooth.vars)==0){
          stop('could not parse smoothing terms s(var,by=phylo) from frmla.phylo')
        } else {
          stop('Current package not suitable for coefficient contrasts of smoothing terms in generalized additive models. To partition smoothing splines, consider algorithm="phylo" or "mStable"')
          warning('Explanatory variables {',paste(smooth.vars,collapse=','),'} defined as partitioning variables for regression coefficient contrasts')
        }
      }
    }
  }

  if (algorithm %in% c('mix','phylo','mStable') & is.null(frmla.phylo)){
    stop('For algorithms "mix", "phylo", and "mStable", must input frmla.phylo')
  }

  if (algorithm %in% c('phylo','mStable')){
    if (!is.null(frmla)){
      warning(paste('input frmla is ignored for algorithm=',algorithm,'. frmla.phylo will be used, frmla be ignored.',sep=''))
    }
    trms <- base::attr(stats::terms(frmla.phylo),'term.labels')
    pvs.in.frmla <- sapply(PartitioningVariables,grepl,trms)
    if (!is.null(dim(pvs.in.frmla))){
      pvs.in.frmla <- apply(pvs.in.frmla,2,any)
    }
    if (any(!pvs.in.frmla)){
      stop(paste('Partitioning Variables {',paste(PartitioningVariables[!pvs.in.frmla],collapse=','),'} not in frmla.phylo',sep=''))
    }
    if (!any(grepl('phylo',trms))){
      stop('term phylo not found as explanatory variable in frmla.phylo')
    }
  }


  output$frmla <- frmla
  output$frmla.phylo <- frmla.phylo

####### Checking Algorithm-input compatibility #######
  if (algorithm=='mStable'){
    mStableAgg=T
  } else {
    mStableAgg=F
  }
  if (algorithm=='mix'){
    if (!(alpha>0 & alpha<=1)){
      stop('alpha must be between 0 and 1')
    }
  }
  if (algorithm!='mStable' & 'matrix' %in% class(Data)){
    stop(paste('algorithm=',algorithm,' requires input "Data" to be of class data.frame or data.table. Input "Data" of class matrix is only allowed for algorithm=mStable',sep=''))
  }
  if ((algorithm=='mStable') & any(grepl('atrix',class(Data))) & is.null(MetaData)){
    if (frmla.phylo[[3]]=='phylo'){
      if (is.null(colnames(Data))){
        MetaData <- data.table('Sample'=1:ncol(Data))
      } else {
        MetaData <- data.table('Sample'=colnames(Data))
      }
    } else {
      stop('algorithm=mStable, with matrix-class input "Data" and frmla.phylo right-hand side containing more variables than just "phylo", requires MetaData input with columns for the RHS of frmla.phylo')
    }
  }
  if ((algorithm=='mStable') & ('data.frame'%in%class(Data)) & (!'Sample' %in% names(Data))){
    stop('data.frame or data.table input "Data" must have column "Sample" for algorithm=mStable')
  }

############## Checking Data and MetaData compatibility ###################
  ## convert Data to data.table
  if (('data.frame' %in% class(Data))&(!'data.table' %in% class(Data))){
    Data <- data.table::as.data.table(Data)
  }

  ## check MetaData for algorithm='mStable'
  ## The main rules are:
  ## (1) if Data is a data.frame or data.table, both MetaData and Data must have "Sample" for alignment.
  ## (2) if Data is a matrix, the columns of Data must match the rows of MetaData.
  ## Otherwise, we'll report the number of samples that overlap.
  if (!is.null(MetaData)){

    if (algorithm!='mStable'){
      warning('input MetaData is ignored for algorithm!="mStable"')
    } else { ## two scenarios: data are matrix/list or data.frame. Either way: check Sample/colnames, make into data.tables & match.
      if (any(class(Data) %in% c('list','matrix'))){ ## matrix input:check/add Sample, data.table MetaData &
        if (class(Data)=='list'){
          if (!all(c(ncol(Data[[1]]),ncol(Data[[2]]))==nrow(MetaData))){
            stop('MetaData does not have same number of rows as columns of Data[[1]] and/or Data[[2]]')
          } else {  #same number of rows/columns - need to match names.
            if (!'Sample' %in% colnames(MetaData)){
              MetaData$Sample <- paste('Sample',1:nrow(MetaData))
            }
            if (is.null(colnames(Data[[1]]))){
              colnames(Data[[1]]) <-  MetaData$Sample
              colnames(Data[[2]]) <- MetaData$Sample
            } else {
              if (!all(colnames(Data[[1]])%in%MetaData$Sample)){
                stop('not all colnames of Data[[1]] are in MetaData$Sample')
              } else {
                Data[[1]] <- Data[[1]][,MetaData$Sample]
                Data[[2]] <- Data[[2]][,MetaData$Sample]
              }
            }
          }
          MetaData <- data.table::as.data.table(MetaData)
          setkey(MetaData,Sample)
          Data[[1]] <- Data[[1]][,MetaData$Sample]
          Data[[2]] <- Data[[2]][,MetaData$Sample]
        } else {
          if (!ncol(Data)==nrow(MetaData)){
            stop('MetaData does not have same number of rows as columns of Data[[1]] and/or Data[[2]]')
          } else {  #same number of rows/columns - need to match names.
            if (!'Sample' %in% colnames(MetaData)){
              MetaData$Sample <- paste('Sample',1:nrow(MetaData))
            }
            if (is.null(colnames(Data))){
              colnames(Data) <-  MetaData$Sample
            } else {
              if (!all(colnames(Data)%in%MetaData$Sample)){
                stop('not all colnames of Data are in MetaData$Sample')
              } else {
                Data <- Data[,MetaData$Sample]
              }
            }
          }
          MetaData <- data.table::as.data.table(MetaData)
          setkey(MetaData,Sample)
          Data <- Data[,MetaData$Sample]
        }

      } else {
        if ((!'Sample' %in% colnames(MetaData)) | (!'Sample' %in% colnames(Data))){
          stop('For input Data of class data.frame and MetaData, both data.frames must have column "Sample" for alignment and aggregation within-samples.')
        }
        if (!'data.table' %in% class(Data)){
          Data <- data.table::as.data.table(Data)
        }
        if (!'data.table' %in% class(MetaData)){
          MetaData <- data.table::as.data.table(MetaData)
        }
        if (!all.equal(sort(unique(MetaData$Sample)),sort(unique(Data$Sample)))){
          warning('Not all Samples are present in both MetaData & Data')
          smpls <- intersect(unique(MetaData$Sample),unique(Data$Sample))
          if (length(smpls)==0){
            stop('No Sample overlap between MetaData and Data')
          } else {
            Data <- Data[Sample %in% smpls]
            MetaData <- MetaData[Sample %in% smpls]
          }
        }
      }
    }
  } else {
    ## meta-data not input. If algorithm='mStable', we must extract MetaData for construction of Data matrix.
    if (algorithm=='mStable'){
      RHS <- setdiff(names(Data),c(as.character(frmla.phylo[[2]]),'phylo','Species'))
      if (!is.null(MetaData)){
        MetaData <- Data[,RHS,with=F]
        MetaData <- data.table:::`[.data.table`(MetaData,!base::duplicated(MetaData))
        data.table::setkey(MetaData,Sample)
      } else {
        if ('list' %in% class(Data)){
          if (is.null(colnames(Data[[1]]))){
            nms <- 1:ncol(Data[[1]])
          } else {
            nms <- colnames(Data[[1]])
          }
        } else {
          if (is.null(colnames(Data))){
            nms <- 1:ncol(Data)
          } else {
            nms <- colnames(Data)
          }
        }
        MetaData <- data.table(Sample=nms)
        
      }
    }
  }

  if ('data.table' %in% class(Data)){
    if (!'Species' %in% names(Data)){
      stop('Data must contain column "Species" whose entries can be found in tip-labels of tree')
    } else {
      if (!all(unique(Data$Species) %in% tree$tip.label)){
        stop('Species in Data not found in tree')
      }
    }
    if (!all(unique(Data$Species) %in% tree$tip.label)){
      stop('Not all species in Data$Species are in tree$tip.label')
    }
    if (!all(tree$tip.label %in% unique(Data$Species))){
      warning('Not all tree$tip.label are in Data$Species - use output pf$tree for downstream analysis')
      tree <- ape::drop.tip(tree,setdiff(tree$tip.label,unique(Data$Species)))
    }
    if (mStableAgg){
      LHS <- setdiff(as.character(frmla.phylo[[2]]),'cbind')
      if (length(LHS)>1){
        if (!all.equal(c('Successes','Failures'),LHS)){
          stop('Unknown left-hand-side of frmla.phylo. For family=binomial, LHS can either be a single variable or cbind(Successes,Failures) - case sensitive.')
        }
      }
      Data <- phyloframe.to.matrix(Data,mat.data = LHS) ##Converts data frame to matrix or list of matrices for each variable on the LHS
      if (length(LHS)>1){
        class(Data) <- c(class(Data),'matrix')
        for (i in 1:length(Data)){
          Data[[i]] <- Data[[i]][tree$tip.label,]
        }
      }
    } else {
      setkey(Data,'Species')
    }
  } else {
    if (class(Data)=='list'){
      if (!base::isTRUE(all.equal(tree$tip.label,rownames(Data[[1]])))){
        if (!all(tree$tip.label %in% rownames(Data[[1]]))){
          stop('Not all tree tip labels are in data')
        } else if (!all(rownames(Data[[1]]) %in% tree$tip.label)){
          stop('Not all rownames of data are in tree tip labels')
        } else {
          Data[[1]] <- Data[[1]][tree$tip.label,]
          Data[[2]] <- Data[[2]][tree$tip.label,]
        }
      }
    } else {
      if (!all(tree$tip.label==rownames(Data))){
        if (!all(tree$tip.label %in% rownames(Data))){
          stop('Not all tree tip labels are in data')
        } else if (!all(rownames(Data) %in% tree$tip.label)){
          stop('Not all rownames of data are in tree tip labels')
        } else {
          Data <- Data[tree$tip.label,]
        }
      }
    }
  }

##### Check additional arguments - family, subset, ... #########
  if ('family' %in% names(list(...))){
    family <- list(...)$family
    if (class(family)=='function'){
      expfamily <- family()$family
    } else {
      expfamily <- family$family
    }
  } else {
    expfamily <- 'gaussian'
  }

  if (algorithm=='mStable'){
    if (expfamily=='binomial'){
      if (!all(c('Successes','Failures') %in% as.character(frmla.phylo[[2]]))){
        stop('Univariate response variables for binomial mStable aggregation not permitted. mStable aggregation for family=binomial requires the left-hand-side of frmla.phylo to have cbind(Successes,Failures).')
      }
      if (!all(c('Successes','Failures') %in% names(Data))){
        if (length(Data)==2 & 'list' %in% class(Data) & all.equal(dim(Data[[1]]),dim(Data[[2]]))){
          warning('Did not name elements of input Data. Assuming the first element is a matrix of Successes and the second element is a matrix of Failures.')
          names(Data) <- c('Successes','Failures')
        } else {
          stop('Unknown input "Data" for family=binomial. Must be either data frame with Successes and Failures, or, for algorithm=mStable, can be a list of two matrices containing Successes and Failures whose rows are species and columns are samples.')
        }
      }
    }
  }

################## Parallelization #############
  if (!is.null(ncores)){
    cl <- phyloFcluster(ncores)
    parallel::clusterExport(cl,varlist = c('Data','MetaData','tree','model.fcn','objective.fcn','cluster.depends'),envir = environment())
    parallel::clusterEvalQ(cl,eval(parse(text=cluster.depends)))
  } else {
    cl <- NULL
  }

############# CoefContrast Calculations ############
  if (algorithm %in% c('CoefContrast','mix')){
    mynorm <- function(b){
      if (length(b)==1){
        return(b^2)
      } else {
        return(matrix(b,nrow=1)%*%matrix(b,ncol=1))
      }
    }
    # species <- unique(Data$Species)
    species <- tree$tip.label
    setkey(Data,Species)
    if ('weights' %in% names(list(...))){
      getModel <- function(spp,Data,frmla,...){
        ix <- Data[Species==spp,which=T]
          return(do.call(model.fcn,args=list('formula'=frmla,'data'=Data,subset=ix,...)))}
    } else {
      getModel <- function(spp,Data,frmla,...){
        return(do.call(model.fcn,args=list('formula'=frmla,'data'=Data[Species==spp],...)))
      }
    }


    tryCatch(models <- lapply(species,getModel,Data,frmla,...),
             error=function(e) stop(paste('Could not implement model.fcn for each species due to following error: \n',e)))
    for (i in 1:length(models)){
      models[[i]]$call <- frmla
    }
    tryCatch(coefficients <- t(sapply(models,stats::coefficients)),
             error=function(e) stop(paste('Could not extract coefficients from model.fcn for each species due to following error \n',e)))
    tryCatch(SE <- t(sapply(models,FUN=function(m) sqrt(diag(stats::vcov(m))))),
             error=function(e) stop(paste('Could not extract standard errors from model.fcn for each species due to following error \n',e)))

    if (all(PartitioningVariables %in% colnames(coefficients))){
      BP <- coefficients[,PartitioningVariables,drop=F]/SE[,PartitioningVariables,drop=F]
    } else {
      pvs_not_found <- PartitioningVariables[!PartitioningVariables %in% colnames(coefficients)]
      pvs_found <- setdiff(PartitioningVariables,pvs_not_found)
      ix <- sapply(pvs_not_found,FUN=function(a,b) grepl(a,b),b=colnames(coefficients)) %>% apply(MARGIN=1,any)
      warning(paste("PartitioningVariables {",paste(pvs_not_found,collapse=','),"} not found in coefficients from glm. Will grep for inexact matches, often caused by PartitioningVariables which are multilevel factors",sep=''))
      warning(paste("Coefficients {",paste(colnames(coefficients)[ix],collapse=','),"} chosen as partitioning variables",sep=''))
      ix <- ix | colnames(coefficients) %in% pvs_found
      if (!any(ix)){
        stop('Could neither match nor grep any PartitioningVariables in the coefficients of model. Try running stats::coefficients on your input model.fcn for a single species to determine the appropriate names for PartitioningVariables.')
      }
      BP <- coefficients[,ix,drop=F]/SE[,ix,drop=F]
    }

    if (nrow(BP)!=length(species)){
      BP <- t(BP)
    }
    rownames(coefficients) <- species
    rownames(BP) <- species
    names(models) <- species
    # models <- models[tree$tip.label]
    # coefficients <- coefficients[tree$tip.label,]
    # BP <- BP[tree$tip.label,]
    output$species.models <- models
    output$coefficient.matrix <- coefficients
    output$coefficient.SE <- SE
    if (algorithm=='mix'){
      rarest.spp <- names(sort(table(Data$Species)))[1]
      ix <- which(tree$tip.label==rarest.spp)
      grp <- list(ix,setdiff(1:length(tree$tip.label),ix))
      tryCatch(getObjective(grp,tree,Data,frmla.phylo,MetaData,PartitioningVariables,mStableAgg=F,expfamily,model.fcn,objective.fcn,min.group.size,...),
               error = function(e) stop(paste('Failure to getObjective from rarest species due to error:',e)))
    }
  }


############# phylofactorization ##############
  treeList <- list(tree)
  binList <- list(1:ape::Ntip(tree))
  Grps=getPhyloGroups(tree)
  pfs=0
  tm <- Sys.time()
  while (pfs<nfactors){

    if (pfs>=1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      Grps <- getNewGroups(tree,treeList,binList)
    }

    if (algorithm %in% c('CoefContrast','mix')){
      V <- t(sapply(Grps,ilrvec,n=length(tree$tip.label)))
      obj <- apply( V %*% BP ,1,mynorm)
      if (algorithm == 'mix'){
        nt <- ceiling(alpha*length(obj))
        Grps <- Grps[order(obj,decreasing = T)[1:nt]]
      }
    }

    if (algorithm!='CoefContrast'){
      if (is.null(ncores)){
        obj <- sapply(Grps,getObjective,
                      tree=tree,
                      Data=Data,
                      frmla=frmla.phylo,
                      MetaData=MetaData,
                      PartitioningVariables=PartitioningVariables,
                      mStableAgg=mStableAgg,
                      expfamily=expfamily,
                      model.fcn=model.fcn,
                      objective.fcn=objective.fcn,
                      min.group.size=min.group.size,
                      ...)
      } else {
        obj <- parallel::parSapply(cl,Grps,
                                   getObjective,
                                   tree=tree,
                                   Data=Data,
                                   frmla=frmla.phylo,
                                   MetaData=MetaData,
                                   PartitioningVariables=PartitioningVariables,
                                   mStableAgg=mStableAgg,
                                   expfamily=expfamily,
                                   model.fcn=model.fcn,
                                   objective.fcn=objective.fcn,
                                   min.group.size=min.group.size,
                                   ...)
      }
    }

    winner <- which.max(obj)
    grp <- Grps[[winner]]
    if (pfs==0){
      if (algorithm=='mStable'){
        output$models <- list(do.call(model.fcn,
                                      args=list('formula'=frmla.phylo,
                                                'data'=mAggregation(Data,grp,tree,MetaData,expfamily,frmla.phylo),...)))
      } else if (algorithm %in% c('phylo','mix')){
        output$models <- list(do.call(model.fcn,args=list('formula'=frmla.phylo,
                                                          'data'=phyloFrame(Data,grp,tree),...)))
      } else {
        output$models <- NULL
      }
    } else {
      if (algorithm=='mStable'){
        output$models <- c(output$models,list(do.call(model.fcn,
                                                      args=list('formula'=frmla.phylo,
                                                                'data'=mAggregation(Data,grp,tree,MetaData,expfamily,frmla.phylo),...))))
      } else if (algorithm %in% c('phylo','mix')){
        output$models <- c(output$models,list(do.call(model.fcn,args=list('formula'=frmla.phylo,
                                                                          'data'=phyloFrame(Data,grp,tree),...))))
      } else {
        output$models <- NULL
      }
    }
    grp <- getLabelledGrp(tree=tree,Groups=Grps[[winner]])
    output$groups <- c(output$groups,list(Grps[[winner]]))

    grpInfo <- matrix(c(names(grp)),nrow=2)
    output$factors <- cbind(output$factors,grpInfo)

    pfs=pfs+1
    tm2 <- Sys.time()
    time.elapsed <- signif(difftime(tm2,tm,units = 'mins'),3)
    if (pfs==1){
      GUI.notification <- paste('\r',pfs,'factor completed in',time.elapsed,'minutes.   ')
    } else {
      GUI.notification <- paste('\r',pfs,'factors completed in',time.elapsed,'minutes.    ')
    }
    if (!is.null(nfactors)){
      GUI.notification <- paste(GUI.notification,'Estimated time of completion:',
                                as.character(tm+difftime(tm2,tm)*nfactors/pfs),
                                '  \r')
    } else {
      GUI.notification <- paste(GUI.notification,'Estimated time of completion: at latest',
                                as.character(tm+difftime(tm2,tm)*nrow(Data)/pfs),
                                '  \r')
    }
    cat(GUI.notification)
    utils::flush.console()
  }

  if (!is.null(ncores)){
    parallel::stopCluster(cl)
    rm('cl')
  }

  if (pfs>1){
    output$factors <- t(output$factors) %>% as.data.frame
  } else {
    output$factors <- matrix(output$factors,nrow=1) %>% as.data.frame
  }
  rownames(output$factors) <- paste('Factor',1:pfs)
  colnames(output$factors) <- c('Group1','Group2')
  output$basis <- matrix(NA,nrow=ape::Ntip(tree),ncol=pfs)
  for (i in 1:length(output$groups)){
    output$basis[,i] <- ilrvec(output$groups[[i]],ape::Ntip(tree))
  }

  output$bins <- bins(output$basis)
  output$tree <- tree
  output$nfactors <- pfs
  if ('data.table' %in% class(Data)){
    if (max(table(Data$Species))==1){
      output$Data <- Data[match(tree$tip.label,Species)]
    } else {
      output$Data <- Data
    }
  } else {
    output$Data <- Data
  }
  if (algorithm!='CoefContrast'){
    for (i in 1:pfs){
      if (!'data.table' %in% class(output$models[[i]])){
        output$models[[i]]$call <- frmla.phylo
      }
    }
  }
  output$MetaData <- MetaData
  output$objective.fcn <- objective.fcn
  output$method <- 'gpf'
  output$PartitioningVariables <- PartitioningVariables
  class(output) <- 'phylofactor'
  output$phylofactor.fcn <- 'gpf'
  return(output)
}
