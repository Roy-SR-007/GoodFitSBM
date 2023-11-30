# Get.Bidirected.Move

# objective :: given a mixed graph G = (d, b); d: directed, b: bidirected, returns an applicable
# bidirected only Markov move in the form of a list	(g.remove, g.add) where g.remove is a graph containing the edges to
# remove and g.add is a graph containing the edges to add; the move may be empty.

# Input::
# d: a directed graph
# b: a bidirected graph

# Output::
# A list of four igraph objects representing the move
# directed igraph object: the directed only edges to remove
# directed igraph object: the directed only edges to add
# undirected igraph object: the reciprocated only edges to remove
# undirected igraph object: the reciprocated only edges to add

#' @import igraph

Get.Bidirected.Move = function(d = NULL, b) {

  if (is.null(d)) {

    d = igraph::graph.empty(vcount(b))

  }

  bidir.piece = Get.Bidirected.Piece(b) # making a call to the `Get.Bidirected.Piece()` routine

  if (is.null(bidir.piece[[1]])) {

    return(list(igraph::graph.empty(vcount(b), directed = FALSE), igraph::graph.empty(vcount(b), directed = FALSE)))

  }
  else {

    # finally, check:
    g.add = igraph::graph(bidir.piece[[2]], n = vcount(b), directed = FALSE)
    g.remove = igraph::graph(bidir.piece[[1]], n = vcount(b),  directed = FALSE)

    # (1) edges to add makes a simple graph; this can happen if more than one partitions
    if (!igraph::is.simple(g.add)) {

      return(list(igraph::graph.empty(vcount(b), directed = FALSE), igraph::graph.empty(vcount(b), directed = FALSE)))

    }

    # (2) edges to add does not intersect b-edges to remove
    # and edges to add does not intersect d-edges to remove
    if (!igraph::ecount(igraph::graph.intersection(igraph::graph.difference(b, g.remove), g.add)) == 0) {

      return(list(igraph::graph.empty(igraph::vcount(b), directed = FALSE), igraph::graph.empty(igraph::vcount(b), directed = FALSE)))

    }

    # (3) neither order of edges to add intersects d
    if (!is.null(d)) {

      if (!igraph::ecount(igraph::graph.intersection(igraph::as.directed(g.add), d)) == 0) {

        return(list(igraph::graph.empty(igraph::vcount(b), directed = FALSE), igraph::graph.empty(igraph::vcount(b), directed = FALSE)))

      }

    }

    return(list(g.remove, g.add))

  }

}
