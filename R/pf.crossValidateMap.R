#' Function for mapping phylogenetic factors to a new community through universal tree
#' 
#' @export
#' @param pf phylofactor class object whose tree tip-labels are all in the \code{universal.tree}
#' @param univeral.tree Universal phylogeny containing all the tip-labels of \code{pf} and all elements of \code{new.community}
#' @param new.community New community
#' @param factors which factors to cross-map
pf.crossValidateMap <- function(pf,universal.tree,new.community,factors=1:pf$nfactors){
  
  
  # universal.tree <- TREE
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
  for (ff in 1:length(factors)){
    which.tree <- which(sapply(treeList,FUN=function(tr,g) all(g %in% tr$tip.label),g=unlist(labelled.grps[[ff]]),simplify=T))
    tr <- treeList[[which.tree]]  #this is the sub-tree being split
    orig.C <- intersect(pf$tree$tip.label,tr$tip.label)
    new.C <- tr$tip.label
    grp <- crossVmap(labelled.grps[[ff]],tr,orig.C,new.C)
    Grps[[ff]] <- lapply(grp,FUN=function(g,new.community) intersect(g,new.community),new.community)
    
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

 return(Grps)
}
