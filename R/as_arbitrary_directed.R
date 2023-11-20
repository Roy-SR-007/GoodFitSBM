# as.arbitrary.directed

# objective :: given an undirected graph b, return a directed graph containing an
# arbitrarily directed copy of each edge of b

# Input::
# b: an undirected graph

# Output::
# a directed graph

as.arbitrary.directed = function(b) {

  # create a directed graph out of the edges of b
  b.decr = igraph::graph(t(igraph::get.edges(b, 1:igraph::ecount(b))))

  # pick a random integer from 0 to no. of edges in b
  num.edges.to.reverse = sample(0:igraph::ecount(b), 1)

  # direct the first num.edges.to.reverse edges in one way and the others the other way
  if (num.edges.to.reverse == 0) {

    b.directed = b.decr

  }
  else {

    random.edge.indices = sample(1:igraph::ecount(b), num.edges.to.reverse)
    b.subset.decr = igraph::graph(t(igraph::get.edges(b, random.edge.indices)))
    el = igraph::get.edgelist(b.subset.decr, names = FALSE) # magically swap cols to reverse direction
    b.subset.incr = igraph::graph(rbind(el[ , 2], el[ , 1]))

    # make the directed graph out of;
    # (reversed.edges.of.b being directed in the decreasing order) union (remaining edges in increasing order)
    b.directed = igraph::graph.union(igraph::graph.difference(b.decr, b.subset.decr), b.subset.incr)

  }

  return(b.directed)

}
