#' Finds interrupting edges along chain of edges
#' @export
#' @param edges non-branching chain of edges
#' @param tree phylo clas object
#' @param original_community original community
#' @param new_community new community - will look for edges from new community interrupting the tree from old community
#' @param old_edges optional input of old edges - the edges of master \code{tree} present in the old community
#' @param new_edges optional input of new edges - the edges of master \code{tree} present in new community.

findInterruptions <- function(edges,tree,original_community,new_community,old_edges=NULL,new_edges=NULL){
  #finds 'interruptions', or sets of taxa which enter into a chain of edges in a tree formed by original_community
  #If repeating this over many sets of edges, can be made more efficient by inputting new_edges and old_edges
  update_new_edges <- FALSE
  if (is.null(old_edges)){
    old_edges <- extractEdges(tree,original_community,type=1)
    update_new_edges <- TRUE
  }
  if (is.null(new_edges) | update_new_edges){
    if (!all(tree$tip.label %in% new_community)){
      new_edges <- extractEdges(tree,new_community,type=1) %>% setdiff(.,old_edges)
    } else {
      new_edges <- setdiff(1:(ape::Nedge(tree)),old_edges)
    }
  }
  
  nds <- c(tree$edge[edges,])
  node.chain <- unique(nds)   #note: these might not be in the right order.
  nds <- as.numeric(names(which(by(nds,nds,length)==1)))
  end.nodes <- as.numeric(names(which(by(nds,nds,length)==1)))
  
  # Only interior nodes can be interrupted
  possible.interruptions <- setdiff(node.chain,end.nodes)
  
  # Let's get a minimal edge matrix for speedy indexing:
  tEdge <- tree$edge
  rownames(tEdge) <- 1:nrow(tEdge)
  
  rt=Ntip(tree)+1
  new_edges_with_root <- new_edges[which(apply(tEdge[new_edges,,drop=F],MARGIN=1,FUN=function(nds,rt) any(nds==rt),rt=rt))]

  ## Now we trim tEdge to contain only the possible interruptions and no other nodes
  tEdge <- apply(tEdge,MARGIN=1,FUN=function(a,possible.interruptions) any(possible.interruptions %in% a),possible.interruptions=possible.interruptions) %>%
              which() %>%
                tEdge[.,]
  
  
  ## Checking possible interruptions and including interrupting groups
  interrupting_grps <- NULL
  interrupting.nds <- NULL
  n <- 0
  for (nd in possible.interruptions){
    #check if node has adjacent edges, not in edge chain and in new_edges
    neighbor_edges <- rownames(tEdge)[which(apply(tEdge,MARGIN=1,FUN=function(a,nd) nd %in% a,nd=nd))] %>% setdiff(.,edges) %>% as.numeric
    interruption <- as.numeric(neighbor_edges[neighbor_edges %in% new_edges])
    if (length(interruption)>0){
      n=n+1
      ## Now, if interruptiong clade does not contain the root of the tree,
      ## the interrupting group is the set of descendants of the neighboring node.
      ## If, however, the interrupting clade contains a new root not present 
      ## in original_community, as in archea being added, then the interrupting group
      ## is the setdiff of new_community and descendants of our interrupting node. 
      if (!any(neighbor_edges %in% new_edges_with_root)){
        # in this case, the descendants of the neighboring node will contain the interrupting group
        neighbor_node <- setdiff(tEdge[toString(interruption),],nd)
        if (neighbor_node<rt){
          interrupting_grps[n] <- list(tree$tip.label[neighbor_node])
        } else {
          interrupting_grps[n] <- list(tree$tip.label[phangorn::Descendants(tree,neighbor_node,type='tips')[[1]]])
        }
      } else { #in this case, the neighbor edges point the way to the root. 
        interrupting_grps[n] <- list(setdiff(tree$tip.label,tree$tip.label[phangorn::Descendants(tree,nd,type='tips')[[1]]]))
      }
      interrupting.nds <- c(interrupting.nds,nd)
    }
    names(interrupting_grps) <- interrupting.nds
  }
  return(interrupting_grps)
}