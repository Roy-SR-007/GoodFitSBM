#'
#' @title Sampling a graph through a Markov move (basis) for ERSBM
#'
#' @description `sample_a_move_ERSBM` to sample a graph in the same fiber; sampling according to the ERSBM (Karwa et al. (2023))
#'
#' @param C a positive integer vector of size n for block assignments of each node; from 1 to K (no of blocks)
#' @param G_current an igraph object which is an undirected graph with no self loop
#'
#' @return A graph
#' \item{sampled graph}{the sampled graph after one move as per the ERSBM}
#'
#' @importFrom igraph graph.empty
#' @importFrom igraph vcount
#' @importFrom igraph graph
#' @importFrom igraph ecount
#' @importFrom igraph graph.intersection
#' @importFrom igraph graph.difference
#' @importFrom igraph as.directed
#' @importFrom igraph is.simple
#' @importFrom igraph is.directed
#' @importFrom igraph graph.union
#' @importFrom igraph get.edges
#' @importFrom igraph get.edge.ids
#' @importFrom igraph as.undirected
#' @importFrom igraph get.edgelist
#' @importFrom igraph subgraph.edges
#' @importFrom igraph E
#' @importFrom igraph V
#' @importFrom igraph graph.complementer
#' @include Get_Next_Network.R
#'
#' @export
#'
#' @seealso [goftest_ERSBM()] performs the goodness-of-fit test for the ERSBM, where graphs are being sampled
#'
#' @examples
#'
#' RNGkind(sample.kind = "Rounding")
#' set.seed(1729)
#'
#' # We model a network with 3 even classes
#' n1 = 50
#' n2 = 50
#' n3 = 50
#'
#' # Generating block assignments for each of the nodes
#' n = n1 + n2 + n3
#' class = rep(c(1, 2, 3), c(n1, n2, n3))
#'
#' # Generating the adjacency matrix of the network
#' # Generate the matrix of connection probabilities
#' cmat = matrix(
#'   c(
#'     30, 0.05, 0.05,
#'     0.05, 30, 0.05,
#'     0.05, 0.05, 30
#'   ),
#'   ncol = 3,
#'   byrow = TRUE
#' )
#' pmat = cmat / n
#'
#' # Creating the n x n adjacency matrix
#' adj <- matrix(0, n, n)
#' for (i in 2:n) {
#'   for (j in 1:(i - 1)) {
#'     p = pmat[class[i], class[j]] # We find the probability of connection with the weights
#'     adj[i, j] = rbinom(1, 1, p) # We include the edge with probability p
#'   }
#' }
#'
#' adjsymm = adj + t(adj)
#'
#' # graph from the adjacency matrix
#' G = igraph::graph_from_adjacency_matrix(adjsymm, mode = "undirected", weighted = NULL)
#'
#' # sampling a Markov move for the ERSBM
#' G_sample = sample_a_move_ERSBM(class, G)
#'
#' # plotting the sampled graph
#' plot(G_sample, main = "The sampled graph after one Markov move for ERSBM")
#'
#' @references
#' Karwa et al. (2023). "Monte Carlo goodness-of-fit tests for degree corrected and related stochastic blockmodels",
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology},
#' <https://doi.org/10.1093/jrsssb/qkad084>


sample_a_move_ERSBM <- function(C, G_current) {
  # Input:
  # G_current: igraph object which is an undirected graph and has no self loop
  # C: numeric vector of size n of block assignment; from 1 to k

  # Getting block information
  num_blocks <- length(unique(C))
  n <- length(C)

  # Getting edge information
  num_edges <- length(igraph::E(G_current))
  all_edges <- igraph::get.edgelist(G_current)
  all_edges <- all_edges[order(all_edges[,1], all_edges[,2]), ]

  # Getting edge information of the complement graph
  G_comp <- igraph::graph.complementer(G_current, loops = FALSE)
  comp_edges <- igraph::get.edgelist(G_comp)
  comp_edges <- comp_edges[order(comp_edges[,1], comp_edges[,2]), ]

  # Determine whether the move is inter-block or intra-block
  type <- sample.int(2, size = 1)

  if (type == 1) {
    # Adding and deleting edges from same block
    # Sample a block
    s <- sample.int(num_blocks, size = 1)

    # Find edges within the sampled block for both the graph and complement graph
    to_delete <- all_edges[((C[as.numeric(all_edges[, 1])] == s) * (C[as.numeric(all_edges[, 2])] == s)) > 0, ]
    to_add <- comp_edges[((C[as.numeric(comp_edges[, 1])] == s) * (C[as.numeric(comp_edges[, 2])] == s)) > 0, ]

    # Getting dimension information
    to_delete_l <- length(to_delete)
    to_add_l <- length(to_add)

    # Check whether the sampled block has at least one edge
    if ((to_delete_l > 0) * (to_add_l > 0)) {
      # Sample an edge to add from complement graph and delete from the graph
      delete_edge <- sample.int(to_delete_l / 2, 1)
      add_edge <- sample.int(to_add_l / 2, 1)

      # If the sampled block has only one edge in the complement graph, then that will be added,
      # otherwise it will add one edge by sampling randomly from complement graph
      if (length(to_add) == 2) {
        G_sample <- igraph::graph.union(G_current, igraph::graph(to_add, n = n, directed = FALSE))
      } else {
        G_sample <- igraph::graph.union(G_current, igraph::graph(to_add[add_edge, ], n = n, directed = FALSE))
      }

      # If the sampled block has only one edge in the graph, then that will be deleted
      # otherwise it will delete one edge by sampling randomly from graph
      if (to_delete_l == 2) {
        G_sample <- igraph::graph.difference(G_sample, igraph::graph(to_delete, n = n, directed = FALSE))
      } else {
        G_sample <- igraph::graph.difference(G_sample, igraph::graph(to_delete[delete_edge, ], n = n, directed = FALSE))
      }

    } else {
      G_sample <- G_current
    }
  } else if (type == 2) {
    # Adding and deleting edges between two different blocks
    # Sample two different blocks
    two_blocks <- sample.int(num_blocks, size = 2, replace = FALSE)
    s <- two_blocks[1]
    t <- two_blocks[2]

    # Find edges between two fixed blocks for both the graph and complement graph
    inter <- all_edges[((C[as.numeric(all_edges[, 1])] == s) * (C[as.numeric(all_edges[, 2])] == t)) + ((C[as.numeric(all_edges[, 1])] == t) * (C[as.numeric(all_edges[, 2])] == s)) > 0, ]
    comp_inter <- comp_edges[((C[as.numeric(comp_edges[, 1])] == s) * (C[as.numeric(comp_edges[, 2])] == t)) + ((C[as.numeric(comp_edges[, 1])] == t) * (C[as.numeric(comp_edges[, 2])] == s)) > 0, ]

    # Getting dimension information
    inter_l <- length(inter)
    comp_inter_l <- length(comp_inter)

    # Check whether the sampled blocks have at least one edge
    if ((inter_l > 0) * (comp_inter_l > 0)) {
      # Sample an edge to add from complement graph and delete from the graph
      delete_edge <- sample.int(inter_l / 2, 1)
      add_edge <- sample.int(comp_inter_l / 2, 1)

      # If the sampled blocks have only one between edge in the complement graph, then that will be added,
      # otherwise it will add one between edge by sampling randomly from complement graph
      if (comp_inter_l == 2) {
        G_sample <- igraph::graph.union(G_current, igraph::graph(comp_inter, n = n, directed = FALSE))
      } else {
        G_sample <- igraph::graph.union(G_current, igraph::graph(comp_inter[add_edge, ], n = n, directed = FALSE))
      }

      # If the sampled blocks have only one between edge in the graph, then that will be added,
      # otherwise it will add one between edge by sampling randomly from graph
      if (inter_l == 2) {
        G_sample <- igraph::graph.difference(G_sample, igraph::graph(inter, n = n, directed = FALSE))
      } else {
        G_sample <- igraph::graph.difference(G_sample, igraph::graph(inter[delete_edge, ], n = n, directed = FALSE))
      }

    } else {
      G_sample <- G_current
    }
  }

  # Output:
  # the graph after one random move
  return(G_sample)
}
