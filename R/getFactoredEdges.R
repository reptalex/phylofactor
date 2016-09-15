#' Returns edges from a phylofactor
#' @export
#' @param tree Phylo object
#' @param grp1 Vector of the indexes for tips contained in group1 (also the nodes corresponding to the tips in group1)
#' @param grp2 same as grp1, but must be a disjoint group from grp1. 
#' @param v Optional input of an ILR-style partition vector, which automatically extracts grp1 and grp2
#' @examples
#' #Case 1: mrca1>root:
#' set.seed(1)
#' tree <- rtree(12)
#' grp1 <- c(9:12)
#' grp2 <- c(1:8) #correct answer: 16
#' getFactoredEdges(tree=tree,grp1=grp1,grp2=grp2)
#' 
#' #Case 2: mrca2>root:
#' grp1 <- c(1:8)
#' grp2 <- c(9:12) #correct answer:16
#' getFactoredEdges(tree=tree,grp1=grp1,grp2=grp2)

#' 
#' grp1 <- c(1:8,12)
#' grp2 <- c(9:11) #correct answer: 17
#' getFactoredEdges(tree=tree,grp1=grp1,grp2=grp2)


#' #Case 3: mrca1>root & mrca2 >root - nodepath passing trough root
#' grp1 <- c(9:11)
#' grp2 <- c(2:3) #correct answer: c(17,16,6,1,3)
#' getFactoredEdges(tree=tree,grp1=grp1,grp2=grp2)


#' #Case 4: mrca1>root & mrca2 >root - nodepath NOT passing trough root
#' grp1 <- c(9:11)
#' grp2 <- c(5,6) #correct answer: c(17,16,7,8,10,11)
#' getFactoredEdges(tree=tree,grp1=grp1,grp2=grp2)


getFactoredEdges <- function(v=NULL,tree,M=NULL,grp1=NULL,grp2=NULL){
  if (is.null(grp1) && is.null(grp2) && is.null(v)){stop('need to input either grp1 and grp2, or ILR vector v')}
  if (is.null(v)==F){
    grp1 <- which(v<0)
    grp2 <- which(v>0)
  }
  
########## the FUNCTION!!
  ntips <- ape::Ntip(tree)
  rooot <- ntips+1

#### get MRCA1
  if (is.null(M)){
    M <- ape::mrca(tree)
  }
if (length(grp1)>1){
  ms <- c(M[grp1,grp1])
  mrca1 <- min(ms[ms>ntips])
} else {
  mrca1 <- grp1
}

#### get MRCA2
if (length(grp2)>1){
  ms <- c(M[grp2,grp2])
  mrca2 <- min(ms[ms>ntips])
  } else {
  mrca2 <- grp2
}
NP <- ape::nodepath(tree,mrca1,mrca2)
#There are a few scenarios: either the groups have mrcas that are both != root
#in which case their shared mrca either (A) is the root or (B) is not the root.
#or, (C) one of the groups has mrca == root. 

# if (A) or (B), we find edges traversing the mrca's of the two groups. 
# if (C), we take the group with MRCA != root, and traverse its nodepath to the root
# until we find the first node containing descendants of the other group.

ix <- which(c(mrca1,mrca2) == rooot)

if (length(ix)!=0){ #only one mrca != root - scenario (C) above. 
  # In this case, the edges separating the two groups are bounded
  # on one side by the mrca of whichever group's mrca<root,
  # and on the other by the first node in NP which has descendants in both grp1 and grp2
  if (ix==2){ #mrca1==root
    # In this case, grp1 we check NP from left to right.
    for (nde in NP[2:length(NP)]){
      descendants <- phangorn::Descendants(tree,nde)[[1]]
      if (any(descendants %in% grp2)){
        NP <- NP[1:which(NP==nde)]
        break
      }
    }
  } else { #mrca2==root
    for (nde in NP[length(NP):2-1]){
      descendants <- phangorn::Descendants(tree,nde)[[1]]
      if (any(descendants %in% grp1)){
        NP <- NP[length(NP):which(NP==nde)]
        break
      }
    }
  }
} else {
  # There's a special case where neither MRCAs are the root, but one clade is contained
  # in the other. We need to check for this case
  rootpath_1 <- ape::nodepath(tree,mrca1,rooot)
  if (mrca2 %in% rootpath_1){
    for (nde in rootpath_1[2:length(rootpath_1)]){
      descendants <- phangorn::Descendants(tree,nde)[[1]]
      if (any(descendants %in% grp2)){
        NP <- rootpath_1[1:which(rootpath_1==nde)]
        break
      }
    }
  } else {
    rootpath_2 <- ape::nodepath(tree,mrca2,rooot)
    if (mrca1 %in% rootpath_2){
      for (nde in rootpath_2[2:length(rootpath_2)]){
        descendants <- phangorn::Descendants(tree,nde)[[1]]
        if (any(descendants %in% grp1)){
          NP <- rootpath_2[1:which(rootpath_2==nde)]
          break
        }
      }
    }
  }
}

L.NP <- length(NP)
edgs <- integer(L.NP-1)
for (nn in 1:(L.NP-1)){
  edgs[nn] <- which(apply(tree$edge,1,function(x,y) all(x %in% y),y=c(NP[nn],NP[nn+1])))
}
return(edgs)
}