# Get.Within.Blocks.Move.beta.SBM

# object :: a bidirected (only) Markov move as the within block move for a beta SBM

# Input, Output:: refer to the `Get.Within.Blocks.Move.beta.SBM()` routine


Get.Within.Blocks.Move.beta.SBM = function(g) {

  return (Get.Bidirected.Move(igraph::graph.empty(vcount(g)), g))

}
