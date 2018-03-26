#' Generalized phylofactorization - currently skeleton for binomial glm; will expand to exponential, gamma etc. soon.
#' 
#' @export
#' @param Data data table containing columns of "Species", 'N' (counts, 0<=N<=3), "Sample". 
#' @param tree phylo class object containing all species in \code{Data}
#' @param X meta-data containing "Sample" and variables found in \code{frmla}
#' @param frmla Formula. If \code{expfamily='binomial'}, must have c(Successes,Failures)~. Otherwise, the variable for data is "Data", e.g. \code{Data~effort}
#' @param PartitioningVariables Character vector containing the variables in \code{frmla} to be used for phylofactorization. Objective function will be the sum of deviance from all variables listed here.
#' @param frmla.phylo Formula used for method "phylo" and "mStable". Can use species-specific effects by multilevel factor "Species", and phylo-factors with "phylo". e.g. cbind(Successes,Failures)~x+Species*z+phylo*y has universal coefficient for x, species-specific coefficients for z, and phylo-factor contrasted coefficients for y.
#' @param nfactors integer for number of factors to find
#' @param ncores integers for number of cores to use for parallelization
#' @param binom.size binom.size of binomial samples for each species. binom.size=1 for presence/absence data.
#' @param expfamily Either "gaussian" or "binomial" - determines the aggregation method.
#' @param model.fcn Regression function, such as glm, gam, glm.nb, gls. Must have column labelled "Deviance" in \code{\link{anova}}.
#' @param objective.fcn Optional input objective function. Takes \code{model.fcn} output as its input, and returns a number which will be maximized to determine phylogenetic factor.
#' @param algorithm Character, either "CoefContrast", "phylo", "mStable" or "mix". "CoefContrast" will partition the standardized coefficient matrix; "phylo" will produce \code{phylo} factors, "mStable" will use \code{phylo} factors for aggregated groups, and "mix" will use coefficient contrasts to identify the top alpha percent of edges and subsequently use the "phylo" algorithm for edge selection.
#' @param alpha Numeric between 0 and 1 (strictly greater than 0), indicating the top quantile of edges to use when \code{algorithm=='mix'}. Default is alpha=0.2
#' @param ... Additional arguments for \code{model.fcn}, e.g. for default \code{\link{glm}}, can use \code{family=binomial} etc.
#' @examples 
#' library(phylofactor)
#' 
#' ilogit <- function(eta) 1/(1+exp(-eta))
#' set.seed(1)
#' m <- 50
#' n <- 200
#' tree <- rtree(m)
#' X <- data.table('y'=rnorm(n),
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
#' eta <- .5*X$z+X$y
#' p <- ilogit(eta)
#' M <- matrix(rbinom(m*n,binom.size,rep(p,times=m)),nrow=m,ncol=n,byrow=T)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- X$Sample
#' 
#' #### the first clade decreases with y ####
#' eta1 <- .5*X$z-X$y
#' p1 <- ilogit(eta1)
#' for (species in clade){
#'    M[species,] <- rbinom(n,binom.size,p1)
#' }
#' #### the second clade weakly decreases with y ####
#' eta2 <- .5*X$z-.3*X$y
#' p2 <- ilogit(eta2)
#' for (species in clade2){
#'    M[species,] <- rbinom(n,binom.size,p2)
#' }
#' 
#' #### Default algorithm: 'mix' #####
#' ### For default can input data matrix or data frame with "Species", "Sample", and all relevant meta-data
#' DF <- matrix.to.phyloframe(M,data.name='Successes')
#' DF[,Failures:=binom.size-Successes]
#' setkey(DF,Sample)
#' DF <- DF[X]
#' 
#' ### DF must have "Species", "Sample", and the LHS of the formula input.
#' ### A separate data frame or data table, X, can have "Sample" and the RHS of the formula.
#' 
#' ### The default algorithm is "mix", which uses CoefContrast to limit the number of edges for selection by algorithm 'phylo'
#' ### This algorithm has the greatest power but is also computationally intensive.
#' ### It's recommended that you input both frmla (used in CoefContrst) and frmla.phylo (used in algorithm 'phylo').
#' ### For species-specific effects in algorithm 'phylo', you can use the variable "Species", e.g.
#' ### frmla.phylo=cbind(Successes,Failures)~Species*z+phylo*y. For universal/shared coefficients for "z" across species, simply use
#' ### frmla.phylo=cbind(Successes,Failures)~z+phylo*y
#' ### Since we're modelling a constant effect of z across species, we want to offset(z) in the formula.
#' ### Let's get the coefficients for that:
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
#' ### Algorithms "phylo", "mix", and "mStable" have a fairly high percent of the computation which is parallelizable.                   
#'
#'                     
#' ### Another algorithm is "CoefContrast". For this algorithm, you need to input the frmla and Partitioning Variables
#'  ### CoefContrast is extremely fast and best done without parallelization (as it is built off matrix multiplication).
#' pf <- gpf(DF,tree,frmla=cbind(Successes,Failures)~z+y,
#'                     PartitioningVariables='y',
#'                     algorithm='CoefContrast',
#'                     family=binomial,
#'                     nfactors=2)
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
#' 
#' ####################### to partition on y, must have phylo* #########
#' pf <- gpf(M,tree,X,frmla.phylo=cbind(Successes,Failures)~z+phylo*y,nfactors=2,
#'           binom.size=binom.size,family=binomial(link='logit'),
#'           PartitioningVariables='y',algorithm='mStable')
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
#' 
#' pf.tree(pf)
#' par(mfrow=c(2,1))
#' phytools::phylo.heatmap(tree,M[,order(X$y)])
#' phytools::phylo.heatmap(tree,ilogit(pf.predict(pf)[,order(X$y)]))
#' 
#' ################# Poisson Regression
#' set.seed(1)
#' eta <- 2*X$z+X$y
#' lambda <- exp(eta)
#' M <- matrix(rpois(m*n,rep(lambda,times=m)),nrow=m,ncol=n,byrow=T)
#' rownames(M) <- tree$tip.label
#' colnames(M) <- X$Sample
#' 
#' #### the first clade decreases with y ####
#' eta1 <- .3*X$z-X$y
#' lambda1 <- exp(eta1)
#' for (species in clade){
#'    M[species,] <- rpois(n,lambda1)
#' }
#' #### the second clade strongly increases with y ####
#' eta2 <- .3*X$z-X$y
#' lambda2 <- exp(eta2)
#' for (species in clade2){
#'    M[species,] <- rpois(n,lambda2)
#' }
#' 
#' 
#' 
#' ##For non-binomial, use "Data" as response variable #########
#' pf <- gpf(M,tree,X,frmla.phylo=Data~phylo*(z+y),nfactors=2,family=poisson,
#'           PartitioningVariables='y',algorithm='mStable')
#' pf$factors
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
#' 
#' par(mfrow=c(2,1))
#' phytools::phylo.heatmap(tree,M[,order(X$y)])
#' phytools::phylo.heatmap(tree,exp(pf.predict(pf)[,order(X$y)]))
#' 
#' ### partition vector of data controlling for sample effort ###
#' set.seed(1)
#' effort <- rnorm(50)
#' eta <- effort-3
#' eta[clade] <- eta[clade]+6
#' eta[clade2] <- eta[clade2]+8
#' Data <- data.table('Species'=tree$tip.label,effort,'Successes'=rbinom(50,1,ilogit(eta)),'Sample'=1)
#' Data$Failures <- 1-Data$Successes
#' pf <- gpf(Data,tree,frmla.phylo=cbind(Successes,Failures)~effort+phylo,
#'           nfactors=2,algorithm='phylo',family=binomial)
#' all.equal(pf$groups[[1]][[1]],clade) & all.equal(pf$groups[[2]][[1]],clade2)
gpf <- function(Data,tree,X=NULL,frmla=NULL,PartitioningVariables=NULL,frmla.phylo=NULL,nfactors=NULL,ncores=NULL,binom.size=1,expfamily='gaussian',model.fcn=stats::glm,objective.fcn=pvDeviance,algorithm='mix',alpha=0.2,...){
  
  output <- NULL
  output$call <- match.call()
  output$model.fcn <- model.fcn
  output$additional.arguments <- list(...)
  output$binom.size=binom.size
  output$algorithm <- algorithm
  output$PartitioningVariables <- PartitioningVariables
  
# Algorithm & frmla checks --------------------------------------------------------
  ########## checking formulas ############
  if (is.null(frmla) & is.null(frmla.phylo)){
    stop('must input either frmla or frmla.phylo')
  }
  if (is.null(frmla) & algorithm %in% c('CoefContrast','mix')){
    stop(paste('must input frmla for algorithm==',algorithm,sep=''))
  }
  
  if (is.null(frmla.phylo)){  
    frmla.terms <- stats::terms(frmla) %>% base::attr('term.labels')
    if (!all(PartitioningVariables %in% c(frmla.terms,''))){
      stop('Some PartitioningVariables are not found in input frmla')
    }
    if (is.null(PartitioningVariables)){
      warning("Did not input PartitioningVariables nor frmla.phylo - default is to partition based on all variables, including (Intercept)")
    }
    if (algorithm == 'mix'){
      warning("Did not input frmla.phylo for algorithm=mix. Default is to assign species-specific coefficients for all but PartitioningVariables")
    }
    if (algorithm %in% c('phylo','mStable','mix')){
      if (is.null(PartitioningVariables)){
        frmla.phylo <- stats::update(frmla,.~phylo*.)
      } else {
        Non_Partitioning_Variables <- setdiff(frmla.terms,PartitioningVariables)
        PartitioningTerms <- sapply(PartitioningVariables,FUN=function(s) paste('phylo*',s,sep='')) %>% paste(collapse='+')
        if (length(unique(Data$Sample))>1){
          Non_Partitioning_Terms <- sapply(Non_Partitioning_Variables,FUN=function(s) paste('Species*',s,sep='')) %>% paste(collapse='+')
        } else {
          Non_Partitioning_Terms <- paste(Non_Partitioning_Variables,collapse='+')
        }
        LHS <- as.character(frmla[[2]])
        if (length(LHS)==1){
          frmla.phylo <- stats::as.formula(paste(paste(LHS,'~',sep=''),paste(Non_Partitioning_Terms,PartitioningTerms,sep='+')))
        } else {
          if (!all.equal(tolower(LHS),c('cbind','successes','failures'))){
            stop('unknown LHS of input formula. Contact alex.d.washburne@gmail.com with error report')
          } else {
            pp <- paste(LHS[1],'(',paste(LHS[2:3],collapse=','),')',sep='')
            frmla.phylo <- stats::as.formula(paste(paste(pp,'~',sep=''),paste(Non_Partitioning_Terms,PartitioningTerms,sep='+'),sep=''))
          }
        }
      }
    }
  }
  output$frmla <- frmla
  output$frmla.phylo <- frmla.phylo
  
  #### Checking Algorithms ####
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
  
  ############## Checking Data and X compatibility ###################
  
  if ((!'data.table' %in% class(Data))&('data.frame' %in% class(Data))){
    Data <- data.table::as.data.table(Data)
  }
  if ('matrix' %in% class(Data) & (!algorithm=='mStable')){
    stop('if algorithm is not "mStable", input "Data" must be data frame or data table containing columns: "Species", "Sample", and the LHS of the input frmla')
  }
  if (!is.null(X)){
    if (!'data.table'%in%class(X)){
      X <- data.table::as.data.table(X)
      if (!'Sample'%in%names(X)){
        stop('Input meta-data X must contain a column named "Sample"')
      }
      setkey(X,Sample)
    }
  } else {
    if (algorithm %in% c('CoefContrast','mix')){
      RHS <- attr(terms(frmla),'term.labels')
    } else {
      RHS <- setdiff(attr(terms(frmla.phylo),'term.labels'),'phylo')
    }
    X <- Data[,c("Sample",RHS),with=F]
    X <- data.table:::`[.data.table`(X,!base::duplicated(X))
    data.table::setkey(X,Sample)
  }
  if ('data.table' %in% class(Data)){
    if (!all(c('Species','Sample') %in% names(Data))){
      stop('Data must contain columns with "Species" and "Sample"')
    } else {
      if (!all(unique(Data$Species) %in% tree$tip.label)){
        stop('Species in Data not found in tree')
      }
    }
    if (!is.null(X)){  
      if (!any(X$Sample %in% Data$Sample)){
        stop('No X$Sample in Data$Sample')
      } else if (!all(X$Sample %in% Data$Sample)){
        warning(paste('Only',length(intersect(X$Sample,Data$Sample)),'samples shared in both Data and X'))
      }
    }
    if (mStableAgg){
      Data <- phylo.frame.to.matrix(Data)
    } else {
      setkey(Data,'Species')
    }
    if (!all(unique(Data$Species) %in% tree$tip.label)){
      stop('Not all species in Data$Species are in tree$tip.label')
    }
    if (!all(tree$tip.label %in% unique(Data$Species))){
      warning('Not all tree$tip.label are in Data$Species - use output pf$tree for downstream analysis')
      tree <- ape::drop.tip(tree,setdiff(tree$tip.label,unique(Data$Species)))
    }
  }
  if ('matrix' %in% class(Data)){
    if (!all(tree$tip.label==rownames(Data))){
      if (!all(tree$tip.label %in% rownames(Data))){
        stop('Not all tree tip labels are in data')
      } else if (!all(rownames(Data) %in% tree$tip.label)){
        stop('Not all rownames of data are in tree tip labels')
      } else {
        Data <- Data[tree$tip.label,]
      }
    }
    if (!is.null(X)){
      if (!any(X$Sample %in% colnames(Data))){
        stop('No X$Sample in colnames(Data)')
      } else if (!all(X$Sample %in% colnames(Data))){
        warning(paste('Only',length(intersect(X$Sample,colnames(Data))),'samples shared in both Data and X'))
      }
    }
  }
  
  if ('family' %in% names(list(...))){
    family <- list(...)$family
    if (class(family)=='function'){
      expfamily <- family()$family
    } else {
      expfamily <- family$family
    }
  }
  
  if (algorithm %in% c('CoefContrast','mix')){
    mynorm <- function(b){
      if (length(b)==1){
        return(b^2)
      } else {
        return(matrix(b,nrow=1)%*%matrix(b,ncol=1))
      }
    }
    species <- unique(Data$Species)
    setkey(Data,Sample)
    Data <- Data[X]
    setkey(Data,Species)
    tryCatch(models <- lapply(species,
                            FUN=function(sp,Data) model.fcn(frmla,data=Data[Species==sp],...),
                            Data=Data), 
             error=function(e) stop(paste('Could not implement model.fcn for each species due to following error: \n',e)))
    tryCatch(coef <- t(sapply(models,stats::coefficients)),
             error=function(e) stop(paste('Could not extract coefficients from model.fcn for each species due to following error \n',e)))
    tryCatch(SE <- t(sapply(models,FUN=function(m) sqrt(diag(stats::vcov(m)))[PartitioningVariables])),
             error=function(e) stop(paste('Could not extrac standard errors from model.fcn for each species due to following error \n',e)))
    
    BP <- coef[,PartitioningVariables]/SE
    if (ncol(BP)!=length(PartitioningVariables)){
      BP <- t(BP)
      if (ncol(BP)!=length(PartitioningVariables)){
        stop('error in converting partitioning variables to coefficient matrix. Perhaps Partitioning variables do not correspond to names of coefficients(model.fcn)')
      }
    }
    rownames(coef) <- species
    rownames(BP) <- species
    names(models) <- species
    models <- models[tree$tip.label]
    coef <- coef[tree$tip.label,]
    BP <- BP[tree$tip.label,]
    output$spp.model.fcns <- models
    output$coefficient.matrix <- coef
    output$coefficient.SE <- SE
    if (algorithm=='mix'){
      # INCOMPLETE
      rarest.spp <- names(sort(table(Data$Species)))[1]
      ix <- which(tree$tip.label==rarest.spp)
      grp <- list(ix,setdiff(1:length(tree$tip.label),ix))
      tryCatch(getObjective(grp,tree,Data,frmla=frmla.phylo,expfamily=expfamily,model.fcn=model.fcn,PartitioningVariables=PartitioningVariables,mStableAgg=F,objective.fcn=objective.fcn,...),
               error = function(e) stop(paste('Failure to getObjective from rarest species due to error:',e)))
    }
  }
  
  if (!is.null(ncores)){
    cl <- phyloFcluster(ncores)
    parallel::clusterEvalQ(cl,library(data.table))
    parallel::clusterExport(cl,varlist = c('Data','X','tree','binom.size','model.fcn'),envir = environment())
  } else {
    cl <- NULL
  }
  
  
  
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
        obj <- sapply(Grps,getObjective,tree,Data,X,binom.size,frmla.phylo,expfamily,model.fcn,PartitioningVariables,mStableAgg,objective.fcn,...)
      } else {
        obj <- parallel::parSapply(cl,Grps,FUN=function(grp,tree,Data,xx,binom.size,frmla.phylo,expfamily,model.fcn,PartitioningVariables,mStableAgg,objective.fcn,...) getObjective(grp,tree,Data,xx,binom.size,frmla.phylo,expfamily,model.fcn,PartitioningVariables,mStableAgg,objective.fcn,...),tree=tree,Data=Data,xx=X,binom.size=binom.size,frmla.phylo=frmla.phylo,expfamily=expfamily,model.fcn=model.fcn,PartitioningVariables=PartitioningVariables,mStableAgg=mStableAgg,objective.fcn=objective.fcn,...)
      }
    }
    
    winner <- which.max(obj)
    grp <- Grps[[winner]]
    if (pfs==0){
      if (algorithm=='mStable'){
        output$models <- list(model.fcn(frmla.phylo,data=mAggregation(Data,grp,tree,X,binom.size,expfamily),...))
      } else if (algorithm %in% c('phylo','mix')){
        output$models <- list(model.fcn(frmla.phylo,data=phyloFrame(Data,grp,tree),...))
      } else {
        output$models <- NA
      }
    } else {
      if (algorithm=='mStable'){
        output$models <- c(output$models,list(model.fcn(frmla.phylo,data=mAggregation(Data,grp,tree,X,binom.size,expfamily),...)))
      } else if (algorithm %in% c('phylo','mix')){
        output$models <- c(output$models,list(model.fcn(frmla.phylo,data=phyloFrame(Data,grp,tree),...)))
      } else {
        output$models <- NA
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
    }
    cat(GUI.notification)
    utils::flush.console()
  }
  
  if (!is.null(ncores)){
    parallel::stopCluster(cl)
    rm('cl')
  }
  
  output$factors <- t(output$factors)
  
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
  output$X <- X
  output$method <- 'gpf'
  class(output) <- 'phylofactor'
  return(output)
}

