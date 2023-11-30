# Get.Within.Blocks.Move.beta.SBM

# objective :: a bidirected (only) Markov move as the within block move for a beta SBM

# Input, Output:: refer to the `Get.Within.Blocks.Move.beta.SBM()` routine

#' @import igraph


Get.Within.Blocks.Move.beta.SBM = function(g) {

  return (Get.Bidirected.Move(igraph::graph.empty(vcount(g)), g))

}
