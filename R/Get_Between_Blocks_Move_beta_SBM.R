# Get.Between.Blocks.Move.beta.SBM

# objective :: extension to the `Get.Directed.Move.p1.ed()' routine corresponding to the beta-SBM

# Input::
# g: an igraph object

# Output::
# return a list of four igraph objects; refer to the 'Get.Directed.Move.p1.ed()' routine

Get.Between.Blocks.Move.beta.SBM = function(g) {

  d = igraph::as.directed(g, mode = c("arbitrary"))
  b = igraph::graph.empty(vcount(d), d = FALSE)

  move = Get.Directed.Move.p1.ed(d, b) # making a call to the routine `Get.Directed.Move.p1.ed()`
  move[[1]] = igraph::as.undirected(move[[1]])
  move[[2]] = igraph::as.undirected(move[[2]])

  return (move)

}
