# Get.Directed.Move.p1.ed

# objective :: given a mixed graph G = (d, b); d: directed and b: bidirected graphs,
# returning a random move removing directed edges from the graph only  applicable to
# the observed network G that is guaranteed to move to a network in the p1 fiber under the
# reciprocation parameter choice specified

# Input::
# d: an igraph (directed graph) object
# b: an igraph (bidirected) object

# Output::
# A list of four igraph objects representing the move,
# directed igraph object: the directed only edges to remove
# directed igraph object: the directed only edges to add
# undirected igraph object: the reciprocated only edges to remove
# undirected igraph object: the reciprocated only edges to add

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
#' @include Get_Directed_Piece.R

Get.Directed.Move.p1.ed = function(d, b) {

  dir.piece = Get.Directed.Piece(d)

  if (is.null(dir.piece[[1]])) {

    return(list(igraph::graph.empty(igraph::vcount(d)), igraph::graph.empty(igraph::vcount(d))))

  }
  else {

    g.remove = igraph::graph(dir.piece[[1]])
    g.add = igraph::graph(dir.piece[[2]])

    # check that edges to add makes a simple graph, and has no bidirected edges
    if (!igraph::is.simple(igraph::as.undirected(g.add, mode = "each"))){

      return(list(igraph::graph.empty(igraph::vcount(d)), igraph::graph.empty(igraph::vcount(d))))

    }

    # check that edges to add does not conflict with d-edges to remove in any direction
    if (!igraph::ecount(igraph::graph.intersection(igraph::as.undirected(igraph::graph.difference(d, g.remove)), igraph::as.undirected(g.add))) == 0){

      return(list(igraph::graph.empty(igraph::vcount(d)), igraph::graph.empty(igraph::vcount(d))))

    }

    # check that edges to add does not conflict with bidirected part of graph
    if (!is.null(b)) {

      if (!igraph::ecount(igraph::graph.intersection(igraph::as.undirected(g.add), b)) == 0) {

        return(list(igraph::graph.empty(igraph::vcount(d)), igraph::graph.empty(igraph::vcount(d))))

      }

    }

    return(list(g.remove, g.add))

  }

}
