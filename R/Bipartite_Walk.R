Bipartite.Walk <- function(edges.to.remove, simple.only=TRUE, multiplicity.bound=NULL) {
  #connect head of (i+1)st edge to tail of ith edge to complete a walk:
  num.edges = nrow(edges.to.remove)
  edges.to.add = c()
  for (i in 1:(num.edges - 1)) {
    edges.to.add = c(edges.to.add, edges.to.remove[i + 1, 1], edges.to.remove[i, 2])
  }
  edges.to.add = c(edges.to.add, edges.to.remove[1, 1], edges.to.remove[num.edges,2])
  if (simple.only){
    # Ensure that edges.to.add form no loops or multiple edges
    if (!is.simple(graph(edges.to.add)))
      return(NULL)
  }
  if (!is.null(multiplicity.bound)){
    #Check that produced edges satisfy given multiplicity bound.
    print("TODO: multiplicity.bound")
    # numvertices = find number of vertices from multiplicity.bound
    # if all(count.multiple(graph( edges.to.add),n=numvertices)>multiplicity.bound)
    # return(NULL)
  }
  return(edges.to.add)
}
