# Bipartite.Walk

# objective :: Given a (randomized) list of edges (edges to remove), it returns a list
# of edges (edges to add) that complete an even closed walk by connecting the endpoints
# of successive edges; this can be thought of as an operation on the parameter graph;
# the simple optional input `simple.only` makes sure only square free moves are produced.

# Input:: a list of edges to be removed

# Output:: returns a list of edges to be added

#' @importFrom igraph graph
#' @importFrom igraph is.simple

Bipartite.Walk = function(edges.to.remove, simple.only = TRUE, multiplicity.bound = NULL) {

  #connect head of (i+1)st edge to tail of ith edge to complete a walk
  num.edges = nrow(edges.to.remove)
  edges.to.add = c()

  for (i in 1:(num.edges - 1)) {

    edges.to.add = c(edges.to.add, edges.to.remove[i + 1, 1], edges.to.remove[i, 2])

  }

  edges.to.add = c(edges.to.add, edges.to.remove[1, 1], edges.to.remove[num.edges, 2])

  if (simple.only) {

    # ensure that edges.to.add form no loops or multiple edges
    if (!igraph::is.simple(igraph::graph(edges.to.add))) {

      return(NULL)

    }

  }

  ## To be done later!
  #if (!is.null(multiplicity.bound)) {

    # check that produced edges satisfy given multiplicity bound

   # print("TODO: multiplicity.bound")
    # numvertices = find number of vertices from multiplicity.bound
    # if all(count.multiple(graph( edges.to.add),n=numvertices)>multiplicity.bound)
    # return(NULL)
  #}
  return(edges.to.add)

}
