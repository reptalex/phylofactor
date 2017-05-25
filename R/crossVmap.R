#' maps groups from old to new communities for phylofactor cross-validation
#' 
#' @export
#' @param Grps list of length two containing the names of taxa (found in tip-labels of \code{tree}) in two groups separated by an edge or chain of edges.
#' @param tree phylo class object. Must contain all \code{Grps}, \code{original_community} and \code{new_community}
#' @param original_community original community being factored (strings corresponding to tip-labels of \code{tree}. This helps us ID the new branches in the \code{new_community}
#' @param new_community new community being cross-validated (strings corresponding to tip-labels of \code{tree})
#' @param ignore.interruptions Logical on whether or not to ignore edges in new-community that interrupt the edges separating \code{Grps}
#' @examples
#' 
#' set.seed(1)
#' tree <- ape::rtree(7)
#'c1 <- c('t1','t2','t5','t6','t3','t4')
#'c2 <- setdiff(tree$tip.label,c('t3','t4','t5'))
#'
#'grp1 <- c('t1','t2')
#'grp2 <- c('t5','t6')
#'Grps <- list(grp1,grp2)
#'
#'factored.edges <- getFactoredEdges(grp1=match(grp1,tree$tip.label),
#'                                   grp2=match(grp2,tree$tip.label),tree=tree)
#'
#'ecols <- rep('black',ape::Nedge(tree))
#'ecols[factored.edges] <- 'red'
#'plot.phylo(tree,use.edge.length = FALSE,type='unrooted',
#'             edge.width = 3,edge.color = ecols)
#'tiplabels(cex=3)
#'# one clade interrupts our chain of branches: t7
#'# another clade: [t3,t4], was factored out previously and is not in our Grps
#'
#'# We have an optional output "ignore.interruptions=T' which omits interruptions
#'# in output of group comparisons. 
#'# if ignore.interruptions=F, we want to output the list of possible groups
#'# for iteration and refinement of which edge along the chain holds an effect.
#'
#'# the default is to ignore interruptions, and it should output
#'# [[t1,t2]],[[t6]]
#'# by succesfully ignoring t7.
#'
#'crossVmap(Grps,tree,c1,c2)
#'
#'#correct
#'
#'#Quick check - let's include two interruptions, one of which contains the root
#'grp1 <- c('t3','t4')
#'grp2 <- c('t5')
#'Grps <- list(grp1,grp2)
#'c1 <- c(grp1,grp2,'t6')
#'c2 <- setdiff(tree$tip.label,'t4')
#'
#'#correct output should be [[t3]],[[t5]]
#'crossVmap(Grps,tree,c1,c2)
#'
#'
#'# Now, suppose we wanted to use those interruptions to refine phylofactorization
#'# to a subset of the original edges
#'
#'crossVmap(Grps,tree,c1,c2,ignore.interruptions = FALSE)
#'
#'#notice the output: the first and last entries are our output from crossVmap
#'#when ignore.interruptions=T. The intervening elements are the interrupting clades
#'#in the order in which they interrupt the factored edges.
#'
crossVmap <- function(Grps,tree,original_community,new_community,ignore.interruptions=T){
  
if (!all(original_community %in% tree$tip.label)){
  stop('Some elements of original_community are not in tree$tip.labels')
}
if (!all(new_community %in% tree$tip.label)){
  stop('Some elements of new_community are not in tree$tip.labels')
}
if (!all(unlist(Grps) %in% original_community)){
  stop('some elements of Grps not in original_community')
}
  
grp1 <- Grps[[1]]
grp2 <- Grps[[2]]

#### get edges
gg1 <- match(grp1,tree$tip.label)
gg2 <- match(grp2,tree$tip.label)
factored.edges <- getFactoredEdges(grp1=gg1,grp2=gg2,tree=tree)



if (length(factored.edges)>1){
  #there's some risk of interruptions
  old_edges <- extractEdges(tree,original_community,type=1)
  if (!all(tree$tip.label %in% new_community)){ #new community does not fill tree
    new_edges <- extractEdges(tree,new_community,type=1) %>% setdiff(.,old_edges)
  } else { #New community does fill tree - this means all edges are in new community
    new_edges <- setdiff(1:(ape::Nedge(tree)),old_edges)
  }
  
  if (length(new_edges)>0){
    interruptions <- findInterruptions(factored.edges,tree,original_community,new_community,old_edges,new_edges)
  } else {
    ignore.interruptions=T
  }
  #If ignore.interruptions==T, we remove all interrupting taxa from our output Grps
  #otherwise, we output groups in order in order (grp1,edge,interruption1,edge,...,edge,grp2)
  #the default setting ignores interruptions
  #by virtue of only considering descendants of the extreme nodes
} else {
  ignore.interruptions=T # no possible interruptions - only one edge in master tree
}

### get terminal nodes for edge-chain:
nds <- c(tree$edge[factored.edges,])
nds <- as.numeric(names(which(by(nds,nds,length)==1)))

### One of three cases will occur: 
#(1) both nodes are on same side of root, and the first node is a descendant of the second
#(2) same, but second node descendant of the first
#(3) nodes are on opposite sides of the tree

# If 1, then grp1 = Descendants of node 1, and grp2 is setdiff
# If 2, vice versa
# If 3, either will work. 

## To ensure our comparisons are in the same direction (e.g. if birds/snakes in our phylofactor, we want birds/snakes in cross-validation)
## We need to ensure grp1 and grp2 line up, e.g. grp1 contains t2, and not grp2 contains t2.

ROOT <- ape::Ntip(tree)+1
np1 <- ape::nodepath(tree,nds[1],ROOT)
if (nds[2] %in% np1){
  g1 <- tree$tip.label[phangorn::Descendants(tree,nds[1],type='tips')[[1]]]
  if (length(factored.edges)>1){
    nd <- np1[match(nds[2],np1)-1]
    g2 <- setdiff(setdiff(tree$tip.label,g1),tree$tip.label[phangorn::Descendants(tree,nd,type='tips')[[1]]])
  } else {
    g2 <- setdiff(tree$tip.label,g1)
  }
} else {
  np2 <- ape::nodepath(tree,nds[2],ROOT)
  if (nds[1] %in% np2){ #nds[1] is the basal node. 
    g1 <- tree$tip.label[phangorn::Descendants(tree,nds[2],type='tips')[[1]]]
    if (length(factored.edges)>1){
      # node before the basal node:
      nd <- np2[match(nds[1],np2)-1]
      g2 <- setdiff(setdiff(tree$tip.label,g1),tree$tip.label[phangorn::Descendants(tree,nd,type='tips')[[1]]])
    } else {
      g2 <- setdiff(tree$tip.label,g1)
    }
    nds <- rev(nds)
  } else {
    # nodes are on opposite side of root
    g1 <- tree$tip.label[phangorn::Descendants(tree,nds[1],type='tips')[[1]]]
    g2 <- tree$tip.label[phangorn::Descendants(tree,nds[2],type='tips')[[1]]]
  }
}

### Quick check: try to align to right grps:
if (!(all(tree$tip.label[gg1] %in% g1))){ #grp1 is not in g1
  if (!all(tree$tip.label[gg2] %in% g1)){
    stop('Unable to align groups for cross-validation')
  } else {
    gg <- g2
    g2 <- g1
    g1 <- gg
    nds <- rev(nds)
  }
}


### Finally: trim the fat:
if (ignore.interruptions){
  grps <- list(g1,g2)
  names(grps) <- c(nds[1],nds[2])
} else {
  ## as it is, g1 starts at nds[1] and g2 starts at nds[2]
  np <- ape::nodepath(tree,nds[1],nds[2])
  ordr <- order(match(names(interruptions),np))
  #sort interrupting grps from g1, I1, I2,...,In, g2
  interruptions <- interruptions[ordr]
  grps <- vector(mode='list',length=length(interruptions)+2)
  grps[[1]] <- g1
  if (length(interruptions)>0){
    grps[2:(length(interruptions)+1)] <- interruptions[1:length(interruptions)]
  }
  grps[[length(interruptions)+2]] <- g2
  names(grps) <- c(nds[1],names(interruptions),nds[2])
}
## We want our grps to correspond to numbers on the 
grps <- lapply(grps,FUN=function(x,new_community) intersect(x,new_community),new_community=new_community)

return(grps)
}