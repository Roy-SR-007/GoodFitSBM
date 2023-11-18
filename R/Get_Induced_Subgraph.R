Get.Induced.Subgraph<-function(g,vertices){
  if (length(vertices)<2)
    return (graph.empty(n=length(vertices), directed=is.directed(g)))
  pairs = combn(vertices,2)
  ei = get.edge.ids(g, pairs)
  ei = ei[ei!=0]
  return (subgraph.edges(g, ei, delete.vertice=FALSE))
}
