# Get.Directed.Piece

#' @import igraph
#' @include Bipartite_Walk.R

Get.Directed.Piece = function(d) {

  # d = directed part of G.
  # pick a random subset E of edges of d and randomly shuffle it
  # (i.e., E = random sample from d of random size)

  if (igraph::ecount(d) == 2) {

    # avoid unwanted behaviour of sample function
    random.subset.of.d = igraph::get.edges(d, 1:2)
    subset.size = 2

  }
  else if (ecount(d) > 2) {

    subset.size = sample(2:igraph::ecount(d), 1) # this is a random integer
    random.edge.indices = sample(1:(igraph::ecount(d)), subset.size)
    random.subset.of.d = igraph::get.edges(d, random.edge.indices)

  }
  else {

    return(NULL)

  }

  # randomly partition E,
  # and for every part E_i, call the routine `Bipartite.Walk(E_i)` at E_i
  # and merge the edges.to.add_i from each of the partitions into a big set edges.to.add

  number.of.partitions = sample(1:(floor(subset.size / 2)), 1)

  # initialize where to store the pieces of the walk:
  edges.to.add = c()
  edges.to.remove = c()
  more.edges = c()
  num.edges.left =subset.size
  s=1 # index

  while(num.edges.left > 1) {

    if (num.edges.left==2){

      k=2 # avoid unwanted behaviour of sample function

    }
    else {

      k = sample(2:num.edges.left, 1) # size of current part

    }

    if (num.edges.left - k == 1){

      k=k+1 # E's assumption on not leaving out that last edge hanging

    }

    more.edges = Bipartite.Walk(random.subset.of.d[s:(s + k - 1), ])

    if (is.null(more.edges)) {

      return(NULL)

    }
    else {

      edges.to.add = c(edges.to.add, more.edges)

    }
    num.edges.left= num.edges.left - k
    s = s + k
  }
  # edges.to.remove has to be in the same format as edges.to.add, so we do the following
  if (!is.null(edges.to.add)) {

    edges.to.remove = as.vector(t(random.subset.of.d))

  }

  return(list(edges.to.remove, edges.to.add))

}
