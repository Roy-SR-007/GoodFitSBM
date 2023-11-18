Get.Between.Blocks.Move.beta.SBM<-function(g){
  d = as.directed(g, mode = c("arbitrary"))
  b = graph.empty(vcount(d), d=FALSE)
  move = Get.Directed.Move.p1.ed(d,b)
  move[[1]] = as.undirected(move[[1]])
  move[[2]] = as.undirected(move[[2]])
  return (move)
}
Get.Within.Blocks.Move.beta.SBM<-function(g){
  return (Get.Bidirected.Move(graph.empty(vcount(g)),g))
}
