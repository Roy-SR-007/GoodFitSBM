# Get.Bidirected.Piece

# objective :: this function computes bidirected move ONLY without checking for conflicts;
# this calls `Bipartite.Walk()` but first checks if edges are a match;
# randomly direct the entire bidirected graph and calls `Get.Directed.Piece()`

#' @import igraph

Get.Bidirected.Piece <- function(b) {

  if (ecount(b) < 2) {

    return(NULL)

  }

  b.directed = as.arbitrary.directed(b)

  return(Get.Directed.Piece(b.directed))

}
