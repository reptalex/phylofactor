
###################### extractEdges ###################################################
###### This function was created by DJ Bennett and downloaded from his GitHub #########
extractEdges <- function(phylo, taxa, type = 1) {
  # Extract edges from a phylo object using 1 of 3 methods
  #
  # Args:
  #  phylo: phylogeny (ape class)
  #  taxa: vector of taxon names
  #  type:
  #     1 -- phylogeny consisting solely of the taxa, default
  #     2 -- edges from taxon tips to terminal node
  #     3 -- edges unique to taxa
  #
  # Return:
  #  vector of edges
  # TODO(01/07/2013): this may be more achievable with a vegan matrix
  if (!type %in% c(1,2,3)) {
    stop("Type must be an integer: 1, 2 or 3.")
  }
  if (!is.vector(taxa) | !is.character(taxa)) {
    stop("Invalid or no taxa given.")
  }
  if (length(taxa) == length (phylo$tip.label)){
    return(phylo$edge)
  }
  if (type == 1 & length (taxa) == 1){
    stop("length(taxa) == 1 :
         Cannot return a single edge for type 1.")
  }
  # start at the tips and step back into the phylogeny ...
  # ...add all connecting edges to a vector...
  # stop when all paths have met at the same node (type = 1)
  # or when all paths have reached the root node (type = 2)
  # or when all the nodes are unique (type = 3)
  edges <- match (match (taxa, phylo$tip.label), phylo$edge[,2])
  end.nodes <- phylo$edge[edges, 1]
  term.node <- length (phylo$tip.label) + 1
  if (all(end.nodes %in% term.node) || (length(taxa)==2 && type==1)) {
    return(edges)
  } else {
    if (type == 3){
      while (any (duplicated (end.nodes))){
        start.node <- end.nodes[duplicated(end.nodes)][1]
        if (sum (phylo$edge[,1] %in% start.node) == sum (end.nodes %in% start.node)){
          edge <- match (start.node, phylo$edge[,2])
          end.node <- phylo$edge[edge,1]
          edges <- c(edges, edge)
          end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        } else {
          end.nodes <- end.nodes[end.nodes != start.node]
        }
      }
    } else {
      while (TRUE){
        end.nodes <- sort (end.nodes, TRUE)
        start.node <- end.nodes[1]
        edge <- match (start.node, phylo$edge[,2])
        end.node <- phylo$edge[edge,1]
        edges <- c(edges, edge)
        end.nodes <- c(end.nodes[!end.nodes %in% start.node], end.node)
        if (type == 2){
          if (sum (term.node == end.nodes) == length (end.nodes)){
            break
          }
        } else {
          if (sum (end.nodes[1] == end.nodes) == length (end.nodes)){
            break
          }
        }
      }
    }
    return (edges)
  }
  }
