#' Performs phylofactorization of a compositional dataset
#' @export
#' @param Data Data matrix whose rows are tip labels of the tree, columns are samples of the same length as X, and whose columns sum to 1
#' @param tree Phylogeny whose tip-labels are row-names in Data.
#' @param X independent variable. If performing multiple regression, X must be a data frame whose columns contain all the independent variables used in \code{frmla}
#' @param frmla Formula for input in GLM. Default formula is Data ~ X. 
#' @param choice Choice, or objective, function for determining the best edges at each iteration. Must be choice='var' or choice='F'. 'var' minimizes residual variance of clr-transformed data, whereas 'F' maximizes the F-statistic from an analysis of variance.
#' @param Grps Optional input of groups to be used in analysis to override the groups used in Tree. for correct format of groups, see output of getGroups
#' @param nfactors Number of clades or factors to produce in phylofactorization. Default, NULL, will iterate phylofactorization until either dim(Data)[1]-1 factors, or until stop.fcn returns T
#' @param quiet Logical, default is \code{FALSE}, indicating whether or not to display standard warnings. 
#' @param trust.me Logical, default \code{FALSE}, indicating whether or not to trust the input Data to be compositional with no zeros.
#' @param stop.fcn Currently, accepts input of 'KS'. Coming soon: input your own function of the environment in phylofactor to determine when to stop.
#' @param stop.early Logical indicating if stop.fcn should be evaluated before (stop.early=T) or after (stop.early=F) choosing an edge maximizing the objective function.
#' @param ncores Number of cores for built-in parallelization of phylofactorization. Parallelizes the extraction of groups, amalgamation of data based on groups, regression, and calculation of objective function. Be warned - this can lead to R taking over a system's memory, so see clusterage for how to control the return of memory from clusters.
#' @param clusterage Age, i.e. number of iterations, for which a phyloFcluster should be used before returning its memory to the system. Default age=Inf - if having troubles with memory, consider a lower clusterage to return memory to system.
#' @param tolerance Tolerance for deviation of column sums of data from 1. if abs(colSums(Data)-1)>tolerance, a warning message will be displayed.
#' @param delta Numerical value for replacement of zeros. Default is 0.65, so zeros will be replaced with 0.65*min(Data[Data>0])
#' @return Phylofactor object, a list containing: "Data", "tree" - inputs from phylofactorization. Output also includes "factors","glms","terminated" - T if stop.fcn terminated factorization, F otherwise - "bins", "bin.sizes", "basis" - basis for projection of data onto phylofactors, and "Monophyletic.Clades" - a list of which bins are monophyletic and have bin.size>1
#' @examples
#'  set.seed(1)
#'  library(ape)
#'  library(phangorn)
#'  library(compositions)
#'
#' tree <- unroot(rtree(20))

#' X <- as.factor(c(rep(0,5),rep(1,5)))
#' sigClades <- Descendants(tree,c(22,28),type='tips')
#' Data <- matrix(rlnorm(20*10,meanlog = 8,sdlog = .5),nrow=20)
#' rownames(Data) <- tree$tip.label
#' colnames(Data) <- X
#' Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
#' Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
#' Data <- t(clo(t(Data)))
#' Bins <- bins(G=sigClades,set=1:20)
#' 
#' phytools::phylo.heatmap(tree,Data)
#' tiplabels()

#' PF <- PhyloFactor(Data,tree,X,nfactors=2)
#' PF$bins
#' all(PF$bins %in% Bins)
#' 
#' #PhyloFactor has built-in parallelization
#' PF.par  <- PhyloFactor(Data,tree,X,nfactors=2,ncores=2)
#' all.equal(PF,PF.par)
#' 
#' #PhyloFactor can also be used for multiple regression
#' 
#' b <- rlnorm(ncol(Data))
#' a <- as.factor(c(rep(0,5),rep(1,5)))
#' X <- data.frame('a'=a,'b'=b)
#' frmla <- Data~a+b
#' PF.M <- PhyloFactor(Data,tree,X,frmla=frmla,nfactors=2)


PhyloFactor <- function(Data,tree,X,frmla = NULL,choice='var',Grps=NULL,nfactors=NULL,quiet=T,trust.me=F,small.output=F,stop.fcn=NULL,stop.early=NULL,KS.Pthreshold=0.01,ncores=NULL,clusterage=Inf,tolerance=1e-10,delta=0.65,smallglm=F,...){
  
  
  #### Housekeeping
  if (typeof(Data) != 'double'){
    warning(' typeof(Data) is not "double" - will coerce with as.matrix(), but recommend using output $data for downstream analysis')
    Data <- as.matrix(Data)
  }
  if (all(rownames(Data) %in% tree$tip.label)==F){stop('some rownames of Data are not found in tree')}
  if (all(tree$tip.label %in% rownames(Data))==F){
    if(!quiet){
      warning('some tips in tree are not found in dataset - output PF$tree will contain a trimmed tree')
    }
    tree <- ape::drop.tip(tree,setdiff(tree$tip.label,rownames(Data)))}
  if (!all(rownames(Data)==tree$tip.label)){
    warning('rows of data are in different order of tree tip-labels - use output$data for downstream analysis, or set Data <- Data[tree$tip.label,]')
    Data <- Data[tree$tip.label,]
  }
  if (is.null(frmla)){frmla=Data ~ X}
  if (choice %in% c('F','var')==F){stop('improper input "choice" - must be either "F" or "var"')}
  if(is.null(nfactors)){nfactors=Inf}
  if(ape::is.rooted(tree)){
    tree <- ape::unroot(tree)}
  
  #### Default treatment of Data ###
  if (!trust.me){
    if (any(Data==0)){
      if (delta==0.65){
        if (!quiet){
          warning('Data has zeros and will receive default modification of zeros. Zeros will be replaced with delta*min(Data[Data>0]), default delta=0.65')
        }
      }
      rplc <- function(x,delta){
        x[x==0]=min(x[x>0])*delta
        return(x)
      }
      
      Data <- apply(Data,MARGIN=2,FUN=rplc,delta=delta)
      
    }
    if (any(abs(colSums(Data)-1)>tolerance)){
      if (!quiet){
        warning('Column Sums of Data are not sufficiently close to 1 - Data will be re-normalized by column sums')
      }
      Data <- t(compositions::clo(t(Data)))
      
      if (any(abs(colSums(Data)-1)>tolerance)){
        if (!quiet){
          warning('Attempt to divide Data by column sums did not bring column sums within "tolerance" of 1 - will proceed with factorization, but such numerical instability may affect the accuracy of the results')
        }
      }
    }
  }
  
  if (small.output){
    warning('For downstream phylofactor functions, you may need to add pf$X, compositional pf$Data <- clo(OTUTable[tree$tip.labels],) and pf$tree (see beginning of PhyloFactor code for how to reproduce these tags). Information on bins can be found by bins(pf$basis).')
  }
  
  treeList <- list(tree)
  binList <- list(1:ape::Ntip(tree))
  nms <- rownames(Data)
  LogData = log(Data)
  
  #### Get list of groups from tree ####
#   if (is.null(Grps)){
#     Grps <- getGroups(tree)
#   }
  
  ix_cl=NULL
  treetips=NULL
  grpsizes=NULL
  tree_map=NULL
  if (is.null(ncores)){
    Grps <- phylofactor::getGroups(tree)
    cl=NULL
  } else {
    
    ############################## Setting up phyloFcluster ################################
    cl <- phyloFcluster(ncores)
    ### To reduce the data transfer to the clusters, we can allocate X and LogData, which remain unchanged throughout.
    ### To pre-allocate memory on the clusters, we can export Y and gg - an ILR vector and glm, respectively.
    Y <- numeric(ncol(Data))
    xx=X
    dataset <- c(list(Y),as.list(xx))
    names(dataset) <- c('Data',names(xx))
    dataset <- model.frame(frmla,data = dataset)
    gg <- glm(frmla,data=dataset)
    parallel::clusterExport(cl,varlist=c('LogData','Y','gg','dataset'),envir=environment())
    
    #### The following variables - treetips, grpsizes, tree_map, ix_cl - change every iteration.
    #### Updated versions will need to be passed to the cluster.
    nms=rownames(Data)
    treetips <- sapply(treeList,FUN=ape::Ntip)
    grpsizes <- sapply(treeList,FUN=function(tree,lg) ape::Nnode(phy=tree,internal.only=lg),lg=F)
    nnodes <- sum(grpsizes)
    cl_node_map <- sample(1:nnodes)  ### By randomizing, we can help clusters have a more even load.
    tree_map <- cumsum(grpsizes) # if tree_map[i-1]<Nde<=tree_map[i], then node is Nde-tree_map[i-1] in tree i.
    ix_cl <- parallel::clusterSplit(cl,cl_node_map)
  }
  
  ### Get OTUs from tree
  OTUs <- tree$tip.label
  if (choice=='var'){
    totalvar= Data %>% apply(.,MARGIN=2,function(x) log(x)-mean(log(x))) %>% apply(.,MARGIN=1,var) %>% sum
  } else {
    totalvar=NULL
  }
  
  ################ OUTPUT ###################
  output <- NULL
  if (!small.output){
    output$Data <- Data
    output$X <- X
    output$tree <- tree
  }
  output$glms <- list()
  n <- length(tree$tip.label)
  
  if (is.null(stop.early) && is.null(stop.fcn)){
    STOP=F
  } else {
    STOP=T
    if (is.null(stop.fcn)){
      stop.fcn='KS'
    }
  }
  
  
  ####### On your marks.... Get set.... #######
  pfs=0
  age=1
  output$terminated=F
  
  
  ########## GO! ##############################
  while (pfs < min(length(OTUs)-1,nfactors)){
    
    
    if (is.null(ncores)==F && age==0){
      cl <- phyloFcluster(ncores)
      parallel::clusterExport(cl,varlist=c('LogData','Y','gg','dataset'),envir=environment())
    }
    
    
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
    PhyloReg <- PhyloRegression(LogData,X,frmla,Grps,choice,treeList,cl,totalvar,ix_cl,treetips,grpsizes,tree_map,quiet,nms,smallglm)
    # PhyloReg <- PhyloRegression(Data,X,frmla,Grps,choice,treeList,cl,totalvar,ix_cl,treetips,grpsizes,tree_map,quiet,nms,...)
    ############################## EARLY STOP #####################################
    ###############################################################################
    
    age=age+1
    
    if (is.null(ncores)==F && age>=clusterage){
      parallel::stopCluster(cl)
      if (exists('cl') && is.null(cl)==F){
        rm(cl)
      }
      gc()
      age=0
    }
    
    if (STOP){
      if (is.null(stop.early)){
        if (!is.null(stop.fcn)){
          if (stop.fcn=='KS'){
            ks <- ks.test(PhyloReg$p.values,'punif')$p.value
            if (ks>KS.Pthreshold){
              output$terminated=T
              break
            }
          } else {
            stop.fcn(PhyloReg)
          }
        }
      }
    }
    
    ############# update output ########################
    if (is.null(ncores)){
      winner <- PhyloReg$winner
      grp <- getLabelledGrp(winner,tree,Grps)
      output$groups <- c(output$groups,list(Grps[[winner]]))
    } else {
      grp <- getLabelledGrp(tree=tree,Groups=PhyloReg$grp,from.parallel=T)
      output$groups <- c(output$groups,list(PhyloReg$grp))
    }
    grpInfo <- matrix(c(names(grp)),nrow=2)
    output$factors <- cbind(output$factors,grpInfo)
    output$glms[[length(output$glms)+1]] <- PhyloReg$glm
    output$basis <- output$basis %>% cbind(PhyloReg$basis)
    if (choice=='var'){
      if (pfs==0){
        output$ExplainedVar <- PhyloReg$explainedvar
      } else {
        output$ExplainedVar <- c(output$ExplainedVar,PhyloReg$explainedvar)
      }
    }
    
    
    ############################## LATE STOP ######################################
    ############# Decide whether or not to stop based on PhyloReg #################
    if (STOP){
      if (!is.null(stop.early)){
        if (!is.null(stop.fcn)){
          if (stop.fcn=='KS'){
            ks <- ks.test(PhyloReg$p.values,'punif')$p.value
            if (ks>KS.Pthreshold){
              output$terminated=T
              break
            }
          } else {
            stop.fcn(.GlobalEnv)
          }
        }
      }
    }
    
    pfs=pfs+1
    gc()
  }
  
  
  
  
  ########### Clean up the Output ####################################
  
  ####### Factors ##########
  if (!is.null(output$factors)){
    pfs <- dim(output$factors)[2]
    colnames(output$factors)=sapply(as.list(1:pfs),FUN=function(a,b) paste(b,a,sep=' '),b='Factor',simplify=T)
    rownames(output$factors)=c('Group1','Group2')
  }
  output$nfactors <- pfs
  output$factors <- t(output$factors)
  pvalues <- sapply(output$glms,FUN=function(gg) summary(aov(gg))[[1]][1,'Pr(>F)'])
  output$factors <- cbind(output$factors,pvalues)
  
  
  
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
  
  

  
  ### ExplainedVar
  if (choice=='var'){
    names(output$ExplainedVar) <- sapply(as.list(1:pfs),FUN=function(a,b) paste(b,a,sep=' '),b='Factor',simplify=T)
  }
  
  
  
  ### Shut down cluster
  if (is.null(ncores)==F && exists('cl')){  #shut down & clean out the cluster before exiting function
    parallel::stopCluster(cl)
    rm(cl)
    gc()
  }
  
  class(output) <- 'phylofactor'
  return(output)
}