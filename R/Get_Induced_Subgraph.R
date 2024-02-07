# Get.Induced.Subgraph

# objective :: method necessary as neither `induced_subgraph()` nor `subgraph.edges()` gives exactly the required result
# alternatively one would need to use named-vertex igraph objects throughout the code

# Input::
# g: an igraph object
# vertices: list of vertices of the graph g

# Output::
# an igraph object representing the subgraph of g induced by vertices and containing all vertices of g

#' @importFrom igraph graph.empty
#' @importFrom igraph vcount
#' @importFrom igraph graph
#' @importFrom igraph ecount
#' @importFrom igraph graph.intersection
#' @importFrom igraph graph.difference
#' @importFrom igraph as.directed
#' @importFrom igraph is.simple
#' @importFrom igraph is.directed
#' @importFrom igraph graph.union
#' @importFrom igraph get.edges
#' @importFrom igraph get.edge.ids
#' @importFrom igraph as.undirected
#' @importFrom igraph get.edgelist
#' @importFrom igraph subgraph.edges
#' @importFrom utils combn

Get.Induced.Subgraph = function(g, vertices) {

  if (length(vertices) < 2) {

    return(igraph::graph.empty(n = length(vertices), directed = igraph::is.directed(g)))
  }

  pairs = utils::combn(vertices, 2)
  ei = igraph::get.edge.ids(g, pairs)
  ei = ei[ei != 0]

  return(igraph::subgraph.edges(g, ei, delete.vertice = FALSE))
}
