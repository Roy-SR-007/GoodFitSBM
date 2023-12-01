# Get.Within.Blocks.Move.beta.SBM

# objective :: a bidirected (only) Markov move as the within block move for a beta SBM

# Input, Output:: refer to the `Get.Within.Blocks.Move.beta.SBM()` routine

#' @importFrom igraph graph.empty
#' @importFrom igraph vcount
#' @include Get_Bidirected_Move.R


Get.Within.Blocks.Move.beta.SBM = function(g) {

  return (Get.Bidirected.Move(igraph::graph.empty(igraph::vcount(g)), g))

}
