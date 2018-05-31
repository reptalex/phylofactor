#' Function for mapping phylogenetic factors in a phylofactor object to a new community through universal tree
#' 
#' @export
#' @param pf phylofactor class object whose tree tip-labels are all in the \code{universal.tree}
#' @param universal.tree Universal phylogeny containing all the tip-labels of \code{pf} and all elements of \code{new.community}
#' @param new.community New community
#' @param factors which factors to cross-map
#' @param fill.empty Logical. Whether or not to fill any empty groups with the species from the original community. If true, the last element of the output will be a data frame indicating the species added to the factor-groups.
#' @examples 
#' library(phylofactor)
#' set.seed(1)
#' m=10
#' n=20
#' tree <- rtree(m)
#' plot(tree,show.tip.label = F)
#' tiplabels()
#' tree$tip.label <- paste('t',1:10,sep='')

#' c1 <- tree$tip.label[setdiff(1:m,c(5,7))]
#' c2 <- tree$tip.label                       #t5 added to group1, t7 ignored
#' c3 <- tree$tip.label[setdiff(1:m,c(4,6))]  #t5 is the only member of our group 1
#' c4 <- tree$tip.label[setdiff(1:m,4:6)]     #our factor will be empty for this community 
#' clade1 <- c(4,6)

#' D <- matrix(0,nrow=m,ncol=n)
#' D[clade1,] <- rep(1:n,each=2)
#' rownames(D) <- tree$tip.label

#' D <- D[c1,]
#' metadata <- data.frame('a'=1:n,'b'=rnorm(n),'c'=rpois(n,2))
#' metadata <- data.frame('a'=1:n,'b'=factor(rep(0,1,each=n/2)))

#' pf <- PhyloFactor(D,tree,metadata,frmla = Data~a,nfactors=1,transform.fcn = I)

#' pf.groupsTospecies(pf)
#' pf.crossValidateMap(pf,tree,c1)  ## same as above
#' pf.crossValidateMap(pf,tree,c2)  ## add t5 to group1
#' pf.crossValidateMap(pf,tree,c3)  ## Group1 contains only 5
#' pf.crossValidateMap(pf,tree,c4)  ## group1 empty
#' 
#' cv <- pf.crossValidateMap(pf,tree,c4,fill.empty = T) 
#' ## group1 same as in our original community - empty groups filled.
#' cv$groups
#' cv$filler.species                                    
#' ## factor1, group1 received species t4, t6 as fillers for an empty group.
pf.crossValidateMap <- function(pf,universal.tree,new.community,factors=1:pf$nfactors,fill.empty=FALSE){
  
  
  universal.community <- unique(c(pf$tree$tip.label,new.community))
  if (!all(universal.community %in% universal.tree$tip.label)){
    stop('Some members of pf$tree$tip.label or new.community are not in universal.tree')
  }
  universal.tree <- drop.tip(universal.tree,setdiff(universal.tree$tip.label,universal.community))
  labelled.grps <- lapply(pf$groups,FUN=function(g,tr) lapply(g,FUN=function(g,tr) tr$tip.label[g],tr=tr),tr=pf$tree)
  
  
  ######## For each factor, we need to iteratively 
  ### (1) find group in universal tree
  ### (2) Convert universal group to a new.community group
  ### (3) Split universal tree - keep track of slipt trees (& bins) with treeList, binList
  ### (4) return Grps
  
  treeList <- list(universal.tree)
  binList <- list(universal.tree$tip.label)
  Grps <- vector(mode='list',length=length(factors))
  filler.species <- NULL
  for (ff in 1:length(factors)){
    which.tree <- which(sapply(treeList,FUN=function(tr,g) all(g %in% tr$tip.label),g=unlist(labelled.grps[[ff]]),simplify=T))
    tr <- treeList[[which.tree]]  #this is the sub-tree being split
    orig.C <- intersect(pf$tree$tip.label,tr$tip.label)
    new.C <- tr$tip.label
    grp <- crossVmap(labelled.grps[[ff]],tr,orig.C,new.C)
    Grps[[ff]] <- lapply(grp,FUN=function(g,new.community) intersect(g,new.community),new.community)
    if (length(Grps[[ff]][[1]])==0 & fill.empty){
      Grps[[ff]][[1]] <- labelled.grps[[ff]][[1]]
      filler.species <- rbind(filler.species,data.frame('factor'=ff,'group'=1,'species'=labelled.grps[[ff]][[1]]))
    }
    if (length(Grps[[ff]][[2]])==0 & fill.empty){
      Grps[[ff]][[2]] <- labelled.grps[[ff]][[2]]
      filler.species <- rbind(filler.species,data.frame('factor'=ff,'group'=2,'species'=labelled.grps[[ff]][[2]]))
    }
    ### update treeList ###############
    tips=sum(sapply(grp,FUN=function(g) length(g)==1)) ## How many tips are there? They will be removed from BinList and recalculated later.
    leaves <- tr$tip.label
    if (tips==0){
      treeList[[which.tree]] <- ape::drop.tip(tr,setdiff(leaves,grp[[1]]))
      treeList[[length(treeList)+1]] <- ape::drop.tip(tr,setdiff(leaves,grp[[2]]))
    } else if (tips==1){
      cld <- which(!sapply(grp,FUN=function(g) length(g)==1))
      treeList[[which.tree]] <- ape::drop.tip(tr,setdiff(leaves,grp[[cld]]))
    } else {
      treeList[[which.tree]] <- NULL
    }
  }

  if (!is.null(filler.species)){
    Grps <- list('groups'=Grps,'filler.species'=filler.species)
  }
 return(Grps)
}
