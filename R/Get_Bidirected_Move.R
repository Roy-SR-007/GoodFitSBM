Get.Bidirected.Move <- function(d=NULL, b) {
  if (is.null(d)){
    d = graph.empty(vcount(b))
  }
  bidir.piece = Get.Bidirected.Piece(b)
  if (is.null(bidir.piece[[1]]))
    return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
  else {
    # Finally, check :
    g.add = graph(bidir.piece[[2]], n=vcount(b), directed = FALSE)
    g.remove = graph(bidir.piece[[1]],n=vcount(b),  directed = FALSE)
    #(1) edges.to.add makes a simple graph. This can happen if more than one partitions.
    if (!is.simple(g.add))
      return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
    #(2) edges.to.add does not intersect b-edges.to.remove [i.e. no conflicts created!]:
    # and edges.to.add does not intersect d-edges.to.remove [i.e. no conflicts created!]:
    if (!ecount(graph.intersection(graph.difference(b, g.remove), g.add)) == 0)
      return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
    #(3) neither order of edges.to.add intersects D:
    if (!is.null(d)) {
      if (!ecount(graph.intersection(as.directed(g.add), d)) == 0)
        return(list(graph.empty(vcount(b), directed=FALSE),graph.empty(vcount(b), directed=FALSE)))
    }
    return(list(g.remove, g.add))
  }
}
