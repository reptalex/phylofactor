#' Performs phylofactorization of a compositional dataset
#'@export
#' @param Data Data matrix whose rows are tip labels of the tree, columns are samples of the same length as X, and whose columns sum to 1
#' @param tree Phylogeny whose tip-labels are row-names in Data. Regardless if phylogeny input is rooted, the output nodes and  will correspond to an unrooted tree.
#' @param X independent variable.
#' @param frmla Formula for input in GLM. Default formula is Data ~ X
#' @param method Method for amalgamating groups and constructing basis. Default method = "ILR" uses the isometric log-ratio transform. Coming soon: method='add' which uses additive grouping with log-ratio regression.
#' @param choice Choice, or objective, function for determining the best edges at each iteration. Must be choice='var' or choice='F'. 'var' minimizes residual variance, whereas 'F' maximizes the F-statistic from an analysis of variance.
#' @param Grps Optional input of groups to be used in analysis to override the groups used in Tree. for correct format of groups, see output of getGroups
#' @param nclades Number of clades or factors to produce in phylofactorization. Default, NULL, will iterate phylofactorization until either dim(Data)[1]-1 factors, or until stop.fcn returns T
#' @param stop.fcn Currently, accepts input of 'KS'. Coming soon: input your own function of the environment in phylofactor to determine when to stop.
#' @param stop.early Logical indicating if stop.fcn should be evaluated before (stop.early=T) or after (stop.early=F) choosing an edge maximizing the objective function.
#' @param ncores Number of cores for built-in parallelization of phylofactorization. Parallelizes the extraction of groups, amalgamation of data based on groups, regression, and calculation of objective function. Be warned - this can lead to R taking over a system's memory, so see clusterage for how to control the return of memory from clusters.
#' @param clusterage Age, i.e. number of iterations, for which a phyloFcluster should be used before returning its memory to the system. Default age=1.
#' @return Phylofactor object, a list containing: "method", "Data", "tree" - inputs from phylofactorization. Output also includes "nodes","glms","terminated" - T if stop.fcn terminated factorization, F otherwise - "atoms", "atom.sizes", "basis" - basis from "method" for projection of data onto phylofactors, and "Monophyletic.Clades" - a list of which atoms are monophyletic and have atom.size>1
#' @examples
#'  set.seed(1)
#' tree <- unroot(rtree(7))

#' X <- as.factor(c(rep(0,5),rep(1,5)))
#' sigClades <- Descendants(tree,c(9,12),type='tips')
#' Data <- matrix(rlnorm(70,meanlog = 8,sdlog = .5),nrow=7)
#' rownames(Data) <- tree$tip.label
#' colnames(Data) <- X
#' Data[sigClades[[1]],X==0] <- Data[sigClades[[1]],X==0]*8
#' Data[sigClades[[2]],X==1] <- Data[sigClades[[2]],X==1]*9
#' Data <- t(clo(t(Data)))

#' PF <- PhyloFactor(Data,tree,X,nclades=2)


PhyloFactor <- function(Data,tree,X,frmla = NULL,method='ILR',choice='var',Grps=NULL,nclades=NULL,stop.fcn=NULL,stop.early=NULL,ncores=NULL,clusterage=1,...){

  #### Housekeeping
  if (all(colnames(Data) != X)){stop('column names of Data do not correspond to X')}
  if (all(rownames(Data) %in% tree$tip.label)==F){stop('some rownames of Data are not found in tree')
     } else {if (all(rownames(Data)!=tree$tip.label)){ #we need to re-arrange the data matrix to correspond to the tree tip labels
       Data <- Data[match(rownames(Data),tree$tip.label),]
     }}
  if (all(tree$tip.label %in% rownames(Data))==F){
    warning('some tips in tree are not found in dataset - output PF$tree will contain a trimmed tree')
    tree <- ape::drop.tip(tree,setdiff(tree$tip.label,rownames(Data)))}
  if (is.null(frmla)){frmla=Data ~ X}
  if (method %in% c('add','ILR')==F){stop('improper input method - must be either "add" or "multiply"')}
  if (choice %in% c('F','var')==F){stop('improper input "choice" - must be either "F" or "var"')}
  if(is.null(nclades)){nclades=Inf}
  if(ape::is.rooted(tree)){
    warning('Tree is rooted. Output nodes will correspond to unrooted, tree <- unroot(tree)')
    tree <- unroot(tree)}


 #### Get list of groups from tree ####
  if(is.null(Grps)){ #inputting groups allows us to make wrappers around PhyloFactor for efficient parallelization.
    Grps <- getGroups(tree)
  }
 # This is a list of 2-element lists containing the partitioning of tips
 # in our tree according to the edges. The groups can be mapped to the tree
 # via the node given names(Grps)[i]. The OTUs corresponding to the Groups can be found with:
 ### Get OTUs from tree
 OTUs <- tree$tip.label

 ################ OUTPUT ###################
 output <- NULL
 output$method <- method
 output$Data <- Data
 output$X <- X
 output$tree <- tree
 n <- length(tree$tip.label)

 if (is.null(stop.early)==F && is.null(stop.fcn)==T){
   stop.fcn='KS'
 }

 pfs=1
 age=0
 output$terminated=F
 cl=NULL

 while (pfs <= min(length(OTUs)-1,nclades)){



   if (pfs>1){
   Data <- t(compositions::clo(t(Data/PhyloReg$residualData)))
   Grps <- removeGroup(Grps,PhyloReg$group) #removeGroup can be made more efficient if need be.
   }


  if (is.null(ncores)==F && age==0){
    cl <- phyloFcluster(ncores)
  }
   ############# Perform Regression on all of Groups, and implement choice function ##############
   # PhyloReg <- PhyloRegression(Data=Data,X=X,frmla=frmla,Grps=Grps,method,choice)
   PhyloReg <- PhyloRegression(Data,X,frmla,Grps,method,choice,cl,...)
   ############################## EARLY STOP #####################################
   ###############################################################################

   age=age+1

   if (is.null(ncores)==F && age>=clusterage){
     stopCluster(cl)
     rm(cl)
     gc()
     age=0
   }

   if (is.null(stop.early)==F){
     if (is.null(stop.fcn)==F){
       if (stop.fcn=='KS'){
         ks <- ks.test(PhyloReg$p.values,runif(length(PhyloReg$p.values)))$p.value
         if (ks>0.01){
           output$terminated=T
           break
           }
       } else {
         stop.fcn(PhyloReg)
       }
     }
   }

   ############# update output ########################
   output$group <- c(output$group,PhyloReg$group)
   output$nodes <- c(output$nodes,PhyloReg$node)
   output$glms <- c(output$glms,PhyloReg$glm)
   if (choice=='var'){
     if (pfs==1){
       output$ExplainedVar <- PhyloReg$explainedvar
     } else {
       output$ExplainedVar <- c(output$ExplainedVar,(1-sum(output$ExplainedVar[1:(pfs-1)]))*PhyloReg$explainedvar)
     }
   }
   output$basis <- output$basis %>% c(PhyloReg$basis) %>% matrix(ncol=pfs,byrow=F)


   ############################## LATE STOP ######################################
   ############# Decide whether or not to stop based on PhyloReg #################
   if (is.null(stop.early)==T){
     if (is.null(stop.fcn)==F){
       if (stop.fcn=='KS'){
         ks <- ks.test(PhyloReg$p.values,runif(length(PhyloReg$p.values)))$p.value
         if (ks>0.01){
           output$terminated=T
           break
         }
       } else {
        stop.fcn(PhyloReg)
       }
     }
   }

   ### If we haven't stopped, let's update the key variables
   pfs=pfs+1
   gc()

 }

 output$atoms <- atoms(output$basis)
 NewOTUs <- output$atoms
 Monophyletic <- unlist(lapply(NewOTUs,function(x,y) return(ape::is.monophyletic(y,x)),y=tree))
 names(output$atoms)[Monophyletic] <- 'Monophyletic'
 names(output$atoms)[!Monophyletic] <- 'Paraphyletic'


 ### Make the atom size distribution data frame ###
 atomsize <- unlist(lapply(NewOTUs,FUN = length))
 ### The atoms are not all OTUs, but vary in size:
 sizes <- as.list(sort(unique(atomsize)))
 nsizes <- unlist(lapply(sizes,FUN=function(x,y){return(sum(y==x))},y=atomsize))
 output$atom.sizes <- data.frame('Atom Size'=unlist(sizes),'Number of Atoms'=nsizes)

 if (choice=='var'){
  names(output$ExplainedVar) <- output$nodes
 }

 if (length(output$nodes)>0){
   if (sum(output$nodes>n)>0){
   names(output$nodes)[output$nodes>n] <- "clade"
   }

   if (sum(output$nodes<=n)>0){
     names(output$nodes)[output$nodes<=n] <- "tip"
   }
 }
 output$Monophyletic.clades <- intersect(which(names(output$atoms)=='Monophyletic'),which(atomsize>1))

 return(output)
}


