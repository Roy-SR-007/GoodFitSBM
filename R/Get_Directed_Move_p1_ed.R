Get.Directed.Move.p1.ed <- function(d,b){
  dir.piece=Get.Directed.Piece(d)
  if (is.null(dir.piece[[1]])){
    return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
  }else{
    g.remove = graph(dir.piece[[1]])
    g.add = graph(dir.piece[[2]])
    # Check that edges.to.add makes a simple graph, and has no bidirected edges.
    if (!is.simple(as.undirected(g.add,mode="each")))
      return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
    # Check that edges.to.add does not conflict with d - edges.to.remove in any direction
    if (!ecount(graph.intersection(as.undirected(graph.difference(d,g.remove)),as.undirected(g.add)))==0)
      return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
    # Check that edges.to.add does not conflict with bidirected part of graph:
    if (!is.null(b)) {
      if (!ecount(graph.intersection(as.undirected(g.add),b))==0){
        return(list(graph.empty(vcount(d)),graph.empty(vcount(d))))
      }
    }
    return (list(g.remove,g.add))
  }
}
