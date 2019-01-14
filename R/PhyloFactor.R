#' Regression-based phylofactorization
#' @export
#' @param Data Data matrix whose rows are tip labels of the tree, columns are samples of the same length as X, and whose columns sum to 1
#' @param tree Phylogeny whose tip-labels are row-names in Data.
#' @param X independent variable. If performing multiple regression, X must be a data frame whose columns contain all the independent variables used in \code{frmla}
#' @param frmla Formula for input in GLM. Default formula is Data ~ X. 
#' @param choice Choice, or objective, function for determining the best edges at each iteration using default regression. Must be choice='var' or choice='F'. 'var' minimizes residual variance of clr-transformed data, whereas 'F' maximizes the F-statistic from an analysis of variance.
#' @param transform.fcn Function for transforming data prior to projection onto contrast bases. Default is \code{log}, in which case zeros are internally replaced by 0.65. The transform function must preserve matrix class objects.
#' @param contrast.fcn Contrast function. Default is an efficient version of \code{BalanceContrast}. Another built-in option is \code{\link{amalgamate}} - for amalgamation-based analyses of compositional data, set \code{transform.fcn=I} and \code{contrast.fcn=amalgamate}.
#' @param method Which default objective function to use either "glm", "max.var" or "gam".
#' @param nfactors Number of clades or factors to produce in phylofactorization. Default, NULL, will iterate phylofactorization until either dim(Data)[1]-1 factors, or until stop.fcn returns T
#' @param small.output Logical, indicating whether or not to trim output. If \code{TRUE}, output may not work with downstream summary and plotting wrappers.
#' @param stop.fcn Currently, accepts input of 'KS'. Coming soon: input your own function of the environment in phylofactor to determine when to stop.
#' @param stop.early Logical indicating if stop.fcn should be evaluated before (stop.early=T) or after (stop.early=F) choosing an edge maximizing the objective function.
#' @param KS.Pthreshold Numeric between 0 and 1. P-value threshold for KS-test as default stopping-function.
#' @param alternative alternative hypothesis input to \code{\link{ks.test}} if KS stopping function is used
#' @param ncores Number of cores for built-in parallelization of phylofactorization. Parallelizes the extraction of groups, amalgamation of data based on groups, regression, and calculation of objective function. Be warned - this can lead to R taking over a system's memory.
#' @param delta Numerical value for replacement of zeros. Default is 0.65, so zeros will be replaced column-wise with 0.65*min(x[x>0])
#' @param choice.fcn Function for customized choice function. Must take as input the numeric vector of ilr coefficients \code{y}, the input meta-data/independent-variable \code{X}, and a logical \code{PF.output}. If \code{PF.output==F}, the output of \code{choice.fcn} must be a two-member list containing numerics \code{output$objective} and \code{output$stopStatistic}. Phylofactor will choose the edge which maximizes \code{output$objective} and a customzed input \code{stop.fcn} can be used with the \code{output$stopStatistics} to stop phylofactor internally. 
#' @param cluster.depends Character parsed and evaluated by cluster to load all dependencies for custom choice.fcn. e.g. \code{cluster.depends <- 'library(bayesm)'}
#' @param ... optional input arguments for \code{\link{glm}} or, if \code{method=='gam'}, input for \code{mgcv::gam}
#' @return Phylofactor object, a list containing: "Data", "tree" - inputs from phylofactorization. Output also includes "factors","glms","terminated" - T if stop.fcn terminated factorization, F otherwise - "bins", "bin.sizes", "basis" - basis for projection of data onto phylofactors, and "Monophyletic.Clades" - a list of which bins are monophyletic and have bin.size>1. For customized \code{choice.fcn}, Phylofactor outputs \code{$custom.output}. 
#' @examples
#' set.seed(2)
#' library(phylofactor)
#' library(phangorn)
#' library(mgcv)
#' mar <- par('mar')
#' clo <- function(X) X/rowSums(X) 
#' 
#' ## Example with pseudo-simulated data: real tree with real taxonomy, but fake abundance patterns.
#' data("FTmicrobiome")
#' tree <- FTmicrobiome$tree
#' Taxonomy <- FTmicrobiome$taxonomy
#' tree <- drop.tip(tree,setdiff(tree$tip.label,sample(tree$tip.label,20)))
#' 
#' ### plot phylogeny ###
#' plot.phylo(tree,use.edge.length=FALSE,main='Community Phylogeny')
#' nodelabels()
#' 
#' Taxonomy <- Taxonomy[match(tree$tip.label,Taxonomy[,1]),]
#' X <- as.factor(c(rep(0,5),rep(1,5)))
#' 
#' ### Simulate data ###
#' Factornodes <- c(37,27)
#' Factoredges <- sapply(Factornodes,FUN=function(n,tree) which(tree$edge[,2]==n),tree=tree)
#' edgelabels(c('PF 1','PF 2'),edge=Factoredges,cex=2,bg='red')
#' sigClades <- Descendants(tree,Factornodes,type='tips')
#' 
#' Data <- matrix(rlnorm(20*10,meanlog = 8,sdlog = .5),nrow=20)
#' rownames(Data) <- tree$tip.label
#' colnames(Data) <- X
#' Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
#' Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
#' Data <- t(clo(t(Data)))
#' Bins <- bins(G=sigClades,set=1:20)
#' pf.heatmap(tree=tree,Data=Data)
#' 
#' ### PhyloFactor ###
#' PF <- PhyloFactor(Data,tree,X,nfactors=2)
#' PF$bins
#' all(PF$bins %in% Bins)
#' 
#' 
#' ######### Summary tools ##########
#' PF$factors                                         
#' # Notice that both of the groups at the first factor are labelled as "Monophyletic"
#' # due to the unrooting of the tree
#' PF$ExplainedVar
#' 
#' # A coarse summary tool
#' s <- pf.summary(PF,Taxonomy,factor=1)   
#' s$group1$IDs                                       # Grabbing group IDs
#' s$group2$IDs
#' 
#' # A tidier summary tool
#' td <- pf.tidy(s)                                 
#' td$`group1, Monophyletic`                          
#' # Simplified group IDs - the unique shortest unique prefixes separating the groups
#' td$`group2, Monophyletic`
#' 
#' ## Plotting with summary tools ##
#' par(mfrow=c(1,1),mar=mar)
#' plot(as.numeric(X),td$`Observed Ratio of group1/group2 geometric means`,
#'        ylab='Average ratio of Group1/Group2',pch=18,cex=2)
#' lines(td$`Predicted ratio of group1/group2`,lwd=2)
#' legend(1,12,legend=c('Observed','Predicted'),pch=c(18,NA),lwd=c(NA,2),
#'        lty=c(NA,1),cex=2)
#' 
#' ######### get and plot Phylogenetic info ####
#' PFedges <- getFactoredEdgesPAR(ncores=2,PF=PF) %>% unlist   
#' ## unlisting is unwise if any factor corresponds to more than one edge
#' PFnodes <- tree$edge[PFedges,2]
#' PFclades <- Descendants(tree,PFnodes,'tips')
#' 
#' par(mfrow=c(3,1))
#' pf.heatmap(tree=tree,Data=Data)
#' # edgelabels(c('Factor 1','Factor 2'),edge=PFedges,bg=c('yellow','red'),cex=2)
#' tiplabels('  ',PFclades[[1]],bg='yellow')
#' tiplabels('  ',PFclades[[2]],bg='red')
#' edgelabels(c('PF1','PF2'),edge=PFedges,bg=c('yellow','red'),cex=2)
#' 
#' ### predicted data matrix given phylofactors
#' pred <- pf.predict(PF)
#' colnames(pred) <- colnames(Data)
#' pf.heatmap(tree=tree,Data=pf.predict(PF))
#' ### residual data
#' resid <- Data/pred
#' resid <- resid %>% t %>% clo %>% t
#' pf.heatmap(tree=tree,Data=resid)
#' par(mar=mar)
#' ##################################
#' 
#' ##################################################
#' ############### Other features: ##################
#' 
#' #### glm-style manipulation of formula, weights, etc. #########
#' #w=1:10
#' #PF.weighted <- PhyloFactor(Data,tree,X,weights=w,nfactors=1)
#' 
#' # predict meta-data with ILR abundances by changing formula & family
#' # PF.predict.X <- PhyloFactor(Data,tree,X,frmla=X~Data,nfactors=2,family=binomial)
#' ### more glm controls: offset, model, subset...
#' #PF.fancy <- PhyloFactor(Data,tree,X,frmla=X~Data,nfactors=2,ncores=2,
#' #family=binomial,weights=w,offset=rnorm(10),model=FALSE,subset=3:8)
#' 
#' #### Stopping Function ###########################
#' PF.stop <- PhyloFactor(Data,tree,X,stop.early=TRUE)
#' PF.stop$terminated 
#' # TRUE - this indicates that the factorization was terminated
#' # when there was sufficiently low signal
#' PF.stop$nfactors   # 2 - the correct number of factors
#' all(PF.stop$bins %in% Bins)   #
#' # TRUE - the factors identified were the correct ones.
#' 
#' #### PhyloFactor has built-in parallelization ####
#' PF.par  <- PhyloFactor(Data,tree,X,nfactors=2,ncores=2)
#' all.equal(PF$factors,PF.par$factors)
#' ##################################################
#' 
#' ######### Phylogenetic PCA - maximize variance ###
#' pf.var <- PhyloFactor(Data,tree,method='max.var',nfactors=2)
#' 
#' ######### Multiple regression ####################
#' b <- rlnorm(ncol(Data))
#' a <- as.factor(c(rep(0,5),rep(1,5)))
#' X <- data.frame('a'=a,'b'=b)
#' frmla <- Data~a+b
#' PF.M <- PhyloFactor(Data,tree,X,frmla=frmla,nfactors=2)
#' PF.M$models[[1]]
#' PF.M.par <- PhyloFactor(Data,tree,X,frmla=frmla,nfactors=2,ncores=2)
#' all.equal(PF.M$factors,PF.M.par$factors)
#' 
#' ####### transform.fcn and contrast.fcn ###########
#' ## If we had Gaussian or approximately Gaussian data, 
#' #GausData <- log(Data)
#' #pf.gaussian <- PhyloFactor(GausData,tree,X,frmla=frmla,nfactors=2,transform.fcn=I)
#' 
#' ## We can also perform amalgamation-style analyses with contrast.fcn
#' #pf.amalg <- PhyloFactor(GausData,tree,X,frmla=frmla,
#' #                        nfactors=2,transform.fcn=I,contrast.fcn=amalgamate)                        
#' ##################################################
#' 
#' 
#' 
#' ############################# CUSTOMIZED CHOICE FUNCTIONS ################################
#' #PhyloFactor can also be used for generalized additive models by inputting choice.fcn 
#' #and cluster.depends to load required packages onto the cluster
#' 
#' ### Let's work with some newly simulated data ####
#' set.seed(1.1)
#' n=100
#' Data <- matrix(rlnorm(20*n,meanlog = 8,sdlog = .5),nrow=20)
#' rownames(Data) <- tree$tip.label
#' a <- rnorm(n)
#' b <- rnorm(n)
#' X <- data.frame(a,b)
#' 
#' Data[sigClades[[1]],] <- t(t(Data[sigClades[[1]],])*(20/(1+exp(5*b)))) 
#' ## This clade has a nonlinear response with b, decreasing for high values of b.
#' 
#' Data[sigClades[[2]],] <- t(t(Data[sigClades[[2]],])*8*a^-2)  
#' ## this clade is abundant only for intermediate values of a.
#' 
#' Data <- t(clo(t(Data)))
#' 
#' par(mfrow=c(2,2))
#' plot(a,gMean(Data[sigClades[[1]],],MARGIN=2),ylab='Group1 gMean')
#' plot(b,gMean(Data[sigClades[[1]],],MARGIN=2),ylab='Group1 gMean')
#' plot(a,gMean(Data[sigClades[[2]],],MARGIN=2),ylab='Group2 gMean')
#' plot(b,gMean(Data[sigClades[[2]],],MARGIN=2),ylab='Group2 gMean')
#' 
#' 
#' ######### To input a custom choice.fcn, it needs to take as input the vector of 
#' ######### ILR coefficients 'y', the input meta-data 'X', and a logical PF.output.
#' ######### The output of the custom choice function when PF.output=T 
#' ######### will be returned in PF$custom.output.
#' 
#' ## Demo choice.fcn - generalized additive modelling ##
#' my_gam <- function(y,X,PF.output=FALSE,...){
#'   dataset <- cbind('Data'=y,X)
#'   gg <- mgcv::gam(Data~s(a)+s(b),data=dataset,...)
#' 
#'   if (PF.output){
#'     return(gg)
#'     break
#'   } else {
#'     output <- NULL
#'     
#'     ## The output of the choice function for PF.output=F must contain two labelled numerics:
#'     ## an "objective" statistic and a "stopStatistics". 
#'     output$objective <- getStats(gg)['ExplainedVar']  
#'     output$stopStatistics <- getStats(gg)['Pval']
#'     return(output)
#'   }
#' }
#' 
#' load.mgcv <- 'library(mgcv)'
#' ######### For parallelization of customized choice function, we may also need to input 
#' ######### cluster.depends which loads all dependencies to cluster.
#' ######### The exact call will be clusterEvalQ(cl,eval(parse(text=cluster.depends)))
#' 
#' 
#' PF.G.par <- PhyloFactor(Data,tree,X,choice.fcn=my_gam,sp=c(1,1),
#'            cluster.depends = load.mgcv,nfactors=2,ncores=2)
#' ######### Or we can use the built-in method='gam' and input e.g. smoothing penalty sp
#' PF.G.par2 <- PhyloFactor(Data,tree,X,method='gam',
#'               frmla=Data~s(a)+s(b),sp=c(1,1),nfactors=2,ncores=2)
#' all(sigClades %in% PF.G.par$bins)
#' PF.G.par$factors
#' 
#' 
#' par(mfrow=c(1,2))
#' for (ff in 1:2){
#'   gm <- PF.G.par$custom.output[[ff]]
#'   grp <- PF.G.par$groups[[ff]]
#'   if (ff==1){
#'    x=b
#'    nd <- X
#'     nd$a <- rep(mean(a),length(a))
#'     pred <- predict(gm,newdata = nd)
#'   } else {
#'     x=a
#'     nd <- X
#'     nd$b <- rep(mean(b),length(b))
#'     pred <- predict(gm,newdata = nd)
#'   }
#'   
#'   y <- BalanceContrast(grp,log(Data))
#'   plot(sort(x),y[order(x)],ylab='ILR Coefficient',
#'         xlab='dominant Independent Variable',
#'         main=paste('Factor',toString(ff),sep=' '))
#'   lines(sort(x),pred[order(x)])
#' }
#' 
#' 
#' ######################## Finding Hutchisonian Niches #####################################
#' ### Example of how to use PhyloFactor to identify Gaussian-shapped Hutchinsonian niches ###
#' set.seed(1)
#' n=1000
#' A <- 20
#' mu=-1
#' sigma=0.9
#' Data <- matrix(rlnorm(20*n,meanlog = 8,sdlog = .5),nrow=20)
#' rownames(Data) <- tree$tip.label
#' X <- rnorm(n)
#' Data[sigClades[[1]],] <- t(t(Data[sigClades[[1]],])*A*exp(-(((X-mu)^2)/(2*sigma^2))))
#' Data <- t(clo(t(Data)))
#' 
#' y1 <- gMean(Data[sigClades[[1]],],MARGIN=2)
#' y2 <- gMean(Data[setdiff(1:20,sigClades[[1]]),],MARGIN=2)
#' ratios <- y1/y2
#' 
#' par(mfrow=c(1,1))
#' plot(X,ratios,
#'  ylab='Group1/Group2 gMean',log='y',
#'  main='Identifying Gaussian-shaped Hutchinsonian Niches',
#'  xlab='Environmental Variable')
#' 
#' frmla=Data~X+I(X^2) 
#' PF.Gaus <- PhyloFactor(Data,tree,frmla=frmla,X,nfactors=1,ncores=2)
#' 
#' all.equal(sigClades[[1]],PF.Gaus$bins[[2]])
#' y <- PF.Gaus$groups[[1]] %>% BalanceContrast(.,log(Data))
#' plot(X,y,ylab='Group1/Group2 gMean',
#'      main='Identifying Gaussian-shaped Hutchinsonian Niches',
#'      xlab='Environmental Variable')
#' lines(sort(X),predict(PF.Gaus$models[[1]])[order(X)],lwd=4,col='green')
#' legend(-2.5,-3,legend=c('Observed','Predicted'),
#'          pch=c(1,NA),col=c('black','green'),lty=c(NA,1),lwd=c(NA,2))
#' 
#' ### Because the regression is performed on an ILR coordinate, getting an estimate 
#' ### about the optimal habitat preference and the width of habitat preferences
#' ### requires a little algebra
#' grp <- PF.Gaus$groups[[1]]
#' r <- length(grp[[1]])
#' s <- length(grp[[2]])
#' coefs <- PF.Gaus$models[[1]]$coefficients
#' a <- coefs['I(X^2)']
#' b <- coefs['X']
#' c <- coefs['(Intercept)']
#' d <- sqrt(r*s/(r+s))
#' sigma.hat <- sqrt(-d/(2*a))
#' mu.hat <- -b/(2*a)
#' A.hat <- exp(c/d+mu.hat^2/(2*sigma.hat^2))
#' names(A.hat) <- NULL
#' names(mu.hat) <- NULL
#' names(sigma.hat) <- NULL
#' c('A'=A,'A.hat'=A.hat)
#' c('mu'=mu,'mu.hat'=mu.hat)             
#' #The optimal environment for this simulated organism is mu=-1
#' c('sigma'=sigma,'sigma.hat'=sigma.hat) #The standard deviation is ~0.9. 

PhyloFactor <- function(Data,tree,X=NULL,frmla = Data~X,choice='var',transform.fcn=log,contrast.fcn=NULL,method='glm',nfactors=NULL,small.output=F,stop.fcn=NULL,stop.early=NULL,KS.Pthreshold=0.01,alternative='greater',ncores=NULL,delta=0.65,choice.fcn=NULL,cluster.depends='',...){
  
  
  ######################################################## Housekeeping #################################################################################
  
  ############ Matching Data to tree
  if (is.null(X)){
    if (!(method=='max.var'|(!is.null(choice.fcn)))){
      stop('Must input X')
    } else {
      if (!is.null(dim(Data))){
        X <- stats::rnorm(ncol(Data))
      } else {
        X <- 1
      }
    }
  }
  if (!all(rownames(Data) %in% tree$tip.label)){stop('some rownames of Data are not found in tree')}
  if (!all(tree$tip.label %in% rownames(Data))){
    warning('some tips in tree are not found in dataset - output PF$tree will contain a trimmed tree')
    tree <- ape::drop.tip(tree,setdiff(tree$tip.label,rownames(Data)))}
  if (!all(rownames(Data)==tree$tip.label)){
    warning('rows of data are in different order of tree tip-labels - use output$data for downstream analysis, or set Data <- Data[output$tree$tip.label,]')
    Data <- Data[tree$tip.label,]
  }
  
  
  ######## checking choice, choice.fcn and cluster.depends
  if (!(choice %in% c('F','var','custom'))){stop('improper input "choice" - must be either "F", "var" or "custom"')}
  if (!is.null(choice.fcn)){
    if (!method=='glm'){warning('Input choice.fcn will override non-default method, i.e. phylofactorization will use choice.fcn.')}
    choice='custom'
    method='glm'
    if (is.null(cluster.depends)){warning('Did not input choice.fcn dependencies - this may cause errors in parallelization due to unavailable dependencies in cluster')}
  } else {
    choice.fcn <- function(y=NULL,X=NULL,PF.output=NULL){
      ch <- NULL
      ch$objective <- 1
      ch$stopStatistics <- 1
      return(ch)
    }
  }
  eval(parse(text=cluster.depends))
  ########### Checking & initializing non-default phylofactorization
  default.phylofactorization= (choice %in% c('F','var') & method=='glm')
  if (!default.phylofactorization){
    ############ Checking method 
    if (method %in% c('max.var','gam')){
      if (method=='max.var'){
        choice.fcn <- phylofactor::VAR
        choice='var'
      } else {
        if (choice=='custom'){
          stop('Cannot use customized choice function for built-in gam. Must be either "F" or "var". Use ? PhyloFactor for help building your own objective function')
        }
        gamchoice=choice
        choice.fcn <- function(y,X,PF.output=FALSE,ff=frmla,gamchoice.=gamchoice,...){
          return(phylofactor::GAM(y,X,PF.output=PF.output,gamfrmla=ff,gamchoice=gamchoice,...))
        }
        cluster.depends <- 'library(mgcv)'
        dataset <- cbind('Data'=numeric(ncol(Data)),X)
        choice='custom'
      }
    }
  }
  
  ###################### Default treatment of Data #################################
  
  if (all.equal(transform.fcn,log)==T){
    if (any(Data)<0){
      stop('For log-transformed data analysis, all entries of Data must be greater than or equal to 0')
    }
    if (any(Data==0)){
      if (delta==0.65){
        warning('Data has zeros and will receive default modification of zeros. Zeros will be replaced column wise with delta*min(x[x>0]), default delta=0.65')
      }
      rplc <- function(x,delta){
        x[x==0]=min(x[x>0])*delta
        return(x)
      }
      Data <- apply(Data,MARGIN=2,FUN=rplc,delta=delta)
    }
  }
  #####################################################################################
  ## small.output warning
  if (small.output){
    warning('For downstream phylofactor functions starting with "pf.", you may need to add pf$X, pf$Data, pf$tree and appropriately-made element "data" to all pf$models.')
  }
  ######################################################## Housekeeping #################################################################################
  
  
  
  ##################################### Pre-Allocation & initialization ##################################
  if(is.null(nfactors)){nfactors=Inf}
  if(ape::is.rooted(tree)){
    tree <- ape::unroot(tree)}
  treeList <- list(tree)
  binList <- list(1:ape::Ntip(tree))
  nms <- rownames(Data)
  
  if (!isTRUE(all.equal(transform.fcn,I))){
    TransformedData = transform.fcn(Data)
  } else {
    TransformedData <- Data
  }
  output <- NULL
  if (!small.output){
    output$Data <- Data
  }
  rm('Data')
  gc()
  
  
  ix_cl=NULL
  treetips=NULL
  grpsizes=NULL
  tree_map=NULL
  if (is.null(ncores)){
    Grps <- phylofactor::getPhyloGroups(tree)
    cl=NULL
  } else {
    
    ############################## Setting up phyloFcluster ################################
    cl <- phyloFcluster(ncores)
    Y <- numeric(ncol(TransformedData))
    gg <- NULL
    if (!(choice == 'custom' | (method %in% c('gam','max.var')))){
        if (is.null(ncol(X))){
          dataset <- c(list(Y),as.list(X))
        } else {
          dataset <- cbind(Y,X)
        }
        names(dataset) <- c('Data',names(X))
        # dataset <- stats::model.frame(frmla,data = dataset)
        if (any(! colnames(X) %in% colnames(dataset))){
          ix <- which(!colnames(X) %in% colnames(dataset))
          for (nn in 1:length(ix)){
          dataset <-cbind(dataset,X[,ix])
          colnames(dataset)[ncol(dataset)] <- colnames(X)[ix[nn]]
          }
        args <- list('data'=dataset,'formula'=frmla,...)
        do.call(stats::glm,args)
      }
    } else {
      ################# export dependencies for choice.fcn ##################################
      dataset=NULL
      parallel::clusterExport(cl,varlist=c('choice.fcn','cluster.depends'),envir=environment())
      parallel::clusterEvalQ(cl,eval(parse(text=cluster.depends)))
      #######################################################################################
    }
    xx <- X
    parallel::clusterExport(cl,varlist=c('xx','X','TransformedData','Y','gg','dataset'),envir=environment())
    
    #### The following variables - treetips, grpsizes, tree_map, ix_cl - change every iteration.
    #### Updated versions will need to be passed to the cluster.
    nms=rownames(TransformedData)
    treetips <- sapply(treeList,FUN=ape::Ntip)
    grpsizes <- sapply(treeList,FUN=function(tree,lg) ape::Nnode(phy=tree,internal.only=lg),lg=F)
    nnodes <- sum(grpsizes)
    cl_node_map <- sample(1:nnodes)  ### By randomizing, we can help clusters have a more even load.
    tree_map <- cumsum(grpsizes) # if tree_map[i-1]<Nde<=tree_map[i], then node is Nde-tree_map[i-1] in tree i.
    ix_cl <- parallel::clusterSplit(cl,cl_node_map)
  }
  
  ### Get OTUs from tree
  OTUs <- tree$tip.label
  n <- length(tree$tip.label)
  if (choice=='var' | method=='max.var'){
    totalvar= TransformedData %>% apply(.,MARGIN=2,function(x) x-mean(x)) %>% apply(.,MARGIN=1,stats::var) %>% sum
  } else {
    totalvar=NULL
  }
  
  ################ OUTPUT Initialization ###################
  if (method=='max.var' | choice=='var'){
    output$total.variance=totalvar
  }
  if (!small.output){
    output$X <- X
    output$tree <- tree
  }
  if (choice != 'custom'){
    if (method != 'max.var'){
      output$models <- list()
    }
  } else {
    output$custom.output <- list()
  }
  output$choice <- choice
  output$choice.fcn <- choice.fcn
  output$cluster.depends <- cluster.depends
  output$method <- method
  output$transform.fcn <- transform.fcn
  output$additional.args <- list(...)
  if (is.null(contrast.fcn)){
    output$contrast.fcn <- BalanceContrast
  } else {
    output$contrast.fcn <- contrast.fcn
  }
  output$stop.fcn <- stop.fcn
  if (is.null(stop.early) && is.null(stop.fcn)){
    STOP=F
  } else {
    STOP=T
    if (is.null(stop.fcn)){
      default.stop=T
    } else {
      default.stop=F
    }
    if (is.null(stop.early)){
      stop.early=T
    } else {
      if (!stop.early){ stop.early=NULL}
    }
  }
  
  
  ####### On your marks.... Get set.... #######
  pfs=0
  output$terminated=F
  tm <- Sys.time()
  ########## GO! ##############################
  while (pfs < min(length(OTUs)-1,nfactors)){
    
    if (pfs>=1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      if (is.null(ncores)){
        Grps <- getNewGroups(tree,treeList,binList)
      } else {
        treetips <- sapply(treeList,FUN=ape::Ntip)
        grpsizes <- sapply(treeList,FUN=function(tree,lg) ape::Nnode(phy=tree,internal.only=lg),lg=F)
        nnodes <- sum(grpsizes)
        cl_node_map <- sample(1:nnodes)  ### By randomizing, we can help clusters have a more even load.
        tree_map <- cumsum(grpsizes) # if tree_map[i-1]<Nde<=tree_map[i], then node is Nde-tree_map[i-1] in tree i.
        ix_cl <- parallel::clusterSplit(cl,cl_node_map)
      }
    }
    
    ############# Perform Regression on all of Groups, and implement choice function ##############
    # PhyloReg <- PhyloRegression(TransformedData,X,frmla,Grps,contrast.fcn,choice,treeList,cl,totalvar,ix_cl,treetips,grpsizes,tree_map,nms,choice.fcn)
    PhyloReg <- PhyloRegression(TransformedData,X,frmla,Grps,contrast.fcn,choice,treeList,cl,totalvar,ix_cl,treetips,grpsizes,tree_map,nms,choice.fcn=choice.fcn,method,...)
    ############################## EARLY STOP #####################################
    ###############################################################################
    
    
    #################### STOP FUNCTIONS ####################
    if (STOP){
      if (!is.null(stop.early)){  #early stop - don't add this factor
          if (default.stop){
            ks <- stats::ks.test(unlist(PhyloReg$stopStatistics),'punif',alternative=alternative)$p.value
            if (ks>KS.Pthreshold){
              if (pfs==0){
                stop('No factors found due to early stop on first iteration')
              }
              output$terminated=T
              break
            }
          } else {
            if(stop.fcn(PhyloReg$stopStatistics)){
              output$terminated=T
              if (pfs==0){
                stop('No factors found due to early stop on first iteration')
              }
              break
            }
          }
      }
    }
    ########################################################
    
    ############# update output ################################################################
    if (is.null(ncores)){
      grp <- getLabelledGrp(tree=tree,Groups=PhyloReg$grp)
      output$groups <- c(output$groups,list(PhyloReg$grp))
    } else {
      grp <- getLabelledGrp(tree=tree,Groups=PhyloReg$grp)
      output$groups <- c(output$groups,list(PhyloReg$grp))
    }
    
    grpInfo <- matrix(c(names(grp)),nrow=2)
    output$factors <- cbind(output$factors,grpInfo)
    
    if (choice != 'custom'){
      if (method != 'max.var'){
        if (small.output){
          PhyloReg$model$data <- NULL
        }
        output$models[[length(output$models)+1]] <- PhyloReg$model
      }
    } else {
      output$custom.output[[length(output$custom.output)+1]] <- PhyloReg$custom.output
    }
    
    output$basis <- output$basis %>% cbind(PhyloReg$basis)
    if (choice=='var'){
      if (pfs==0){
        output$ExplainedVar <- PhyloReg$explainedvar
      } else {
        output$ExplainedVar <- c(output$ExplainedVar,PhyloReg$explainedvar)
      }
    }
    ###########################################################################################
    
    ############################## LATE STOP ###############
    #################### STOP FUNCTIONS ####################
    if (STOP){
      if (is.null(stop.early)){
          if (default.stop){
            ks <- stats::ks.test(unlist(PhyloReg$stopStatistics),'punif',alternative=alternative)$p.value
            if (ks>KS.Pthreshold){
              output$terminated=T
              break
            }
          } else {
            if(stop.fcn(PhyloReg$stopStatistics)){
              output$terminated=T
              break
            }
          }
      }
    }
    ########################################################
    
    pfs=pfs+1
    gc()
    tm2 <- Sys.time()
    time.elapsed <- signif(difftime(tm2,tm,units = 'mins'),3)
    if (pfs==1){
      GUI.notification <- paste('\r',pfs,'factor completed in',time.elapsed,'minutes.   ')
    } else {
      GUI.notification <- paste('\r',pfs,'factors completed in',time.elapsed,'minutes.    ')
    }
    if (!is.infinite(nfactors)){
      GUI.notification <- paste(GUI.notification,'Estimated time of completion:',
                                as.character(tm+difftime(tm2,tm)*nfactors/pfs),
                                '  \r')
    } else {
      GUI.notification <- paste(GUI.notification,'Estimated time of completion: at latest',
                                as.character(tm+difftime(tm2,tm)*nrow(TransformedData)/pfs),
                                '  \r')
    }
    base::cat(GUI.notification)
    utils::flush.console()
  }
  
  
  
  
  ########### Clean up the Output ####################################
  
  ####### Factors ##########
  if (!is.null(output$factors)){
    pfs <- dim(output$factors)[2]
    colnames(output$factors)=sapply(as.list(1:pfs),FUN=function(a,b) paste(b,a,sep=' '),b='Factor',simplify=T)
    rownames(output$factors)=c('Group1','Group2')
  }
  output$nfactors <- pfs
  output$factors <- t(output$factors) %>% as.data.frame
  
  if (choice != 'custom'){
    if (method != 'max.var'){
      summary.statistics <- sapply(output$models,FUN=function(gg) getStats(gg)) %>% t %>% as.data.frame()
      colnames(summary.statistics) <- c('Pr(>F)','F','ExpVar')
      summary.statistics <- summary.statistics[,c(3,2,1)]
    } else {
      summary.statistics <- data.frame('ExpVar'=numeric(output$nfactors))
    }
    summary.statistics$ExpVar <- output$ExplainedVar
    output$factors <- cbind(output$factors,summary.statistics)
  }
  
  
  
  ####### Bins ###########
  if (!small.output){
    output$bins <- bins(output$basis)
    NewOTUs <- output$bins
    Monophyletic <- unlist(lapply(NewOTUs,function(x,y) return(ape::is.monophyletic(y,x)),y=tree))
    names(output$bins)[Monophyletic] <- 'Monophyletic'
    names(output$bins)[!Monophyletic] <- 'Paraphyletic'
    ### Make the bin size distribution data frame ###
    binsize <- unlist(lapply(NewOTUs,FUN = length))
    ### The bins are not all OTUs, but vary in size:
    sizes <- as.list(sort(unique(binsize)))
    nsizes <- unlist(lapply(sizes,FUN=function(x,y){return(sum(y==x))},y=binsize))
    output$bin.sizes <- data.frame('Bin Size'=unlist(sizes),'Number of Bins'=nsizes)
    output$Monophyletic.clades <- intersect(which(names(output$bins)=='Monophyletic'),which(binsize>1))
  }
  
  ### Groups
  
  names(output$groups) <- sapply(as.list(1:length(output$groups)),FUN=function(a) paste('factor',toString(a)))
  for (nn in 1:length(output$groups)){
    names(output$groups[[nn]]) <- c('Group1','Group2')
  }
  
  if (output$nfactors>1){
    rownames(output$basis) <- tree$tip.label
  } else {
    names(output$basis) <- tree$tip.label
  }

  
  ### ExplainedVar
  if (choice=='var'){
    names(output$ExplainedVar) <- sapply(as.list(1:pfs),FUN=function(a,b) paste(b,a,sep=' '),b='Factor',simplify=T)
  }
  if (method=='max.var'){
    names(output)[names(output)=='custom.output']='ExplainedVar'
    output$ExplainedVar <- unlist(output$ExplainedVar)/output$total.variance
  }
  if (method=='max.var'){
    output$X <- NULL
  }
  if (method=='gam'){
    output$choice <- gamchoice
    output$models <- output$custom.output
    output$cluster.depends <- NULL
    output$choice.fcn <- NULL
  }
  
  if (length(output$models)>0){
    for (i in 1:output$nfactors){
      output$models[[i]]$call <- frmla
    }
  }
  
  ### Shut down cluster
  if (is.null(ncores)==F && exists('cl')){  #shut down & clean out the cluster before exiting function
    parallel::stopCluster(cl)
    rm(cl)
    gc()
  }
  
  output$phylofactor.fcn <- 'PhyloFactor'
  class(output) <- 'phylofactor'
  return(output)
}
