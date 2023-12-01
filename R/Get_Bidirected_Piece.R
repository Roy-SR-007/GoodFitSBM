# Get.Bidirected.Piece

# objective :: this function computes bidirected move ONLY without checking for conflicts;
# this calls `Bipartite.Walk()` but first checks if edges are a match;
# randomly direct the entire bidirected graph and calls `Get.Directed.Piece()`

#' @importFrom igraph graph.empty
#' @importFrom igraph vcount
#' @importFrom igraph graph
#' @importFrom igraph ecount
#' @importFrom igraph graph.intersection
#' @importFrom igraph graph.difference
#' @importFrom igraph as.directed
#' @importFrom igraph is.simple
#' @importFrom igraph graph.union
#' @importFrom igraph get.edges
#' @importFrom igraph as.undirected
#' @importFrom igraph get.edgelist
#' @include as_arbitrary_directed.R
#' @include Get_Directed_Piece.R

Get.Bidirected.Piece <- function(b) {

  if (igraph::ecount(b) < 2) {

    return(NULL)

  }

  b.directed = as.arbitrary.directed(b)

  return(Get.Directed.Piece(b.directed))

}
