#'
#' @title Maximum Likelihood Estimation of edge probabilities between blocks of a graph
#'
#' @description `get_mle_ERSBM` obtains MLE for the probability of edges between blocks in a graph, used in calculating the goodness-of-fit test statistic for the ERSBM (Karwa et al. (2023))
#'
#' @param G an igraph object which is an undirected graph with no self loop
#' @param C a positive integer vector of size n for block assignments of each node; from 1 to K (no of blocks)
#'
#' @return A matrix of maximum likelihood estimates
#' \item{mleMatr}{a matrix containing the estimated edge probabilities between blocks in a graph}
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
#' @importFrom stats loglin
#'
#' @export
#'
#' @seealso [goftest()] performs the goodness-of-fit test for the beta SBM, where the MLE of the edge probabilities are required
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
#' # mle of the edge probabilities
#' get_mle_ERSBM(G, class)
#'
#' @references
#' Karwa et al. (2023). "Monte Carlo goodness-of-fit tests for degree corrected and related stochastic blockmodels",
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology},
#' <https://doi.org/10.1093/jrsssb/qkad084>

get_mle_ERSBM <- function(G, C) {
  # Input:
  # G_current: G_obs igraph object which is an undirected graph and has no self loop
  # C: numeric vector of size n of block assignment; from 1 to k

  # Getting graph information
  k <- length(unique(C)) # no of block
  n <- length(C) # no of nodes
  A <- igraph::get.adjacency(G, type = "both")

  # Calculate total observed edges between block and within block
  table_obs <- matrix(0, nrow = k, ncol = k)
  for (i in 1:k) {
    for (j in 1:k) {
      table_obs[i, j] <- sum(A[C == i, C == j])
    }
  }

  # Calculate total possible edges between and within block based on nodes information
  tot_edge <- tcrossprod(as.vector(table(C))) - diag(as.vector(table(C)))

  # Estimated value for q_i,j in matrix form
  mle_mat <- table_obs / tot_edge

  # Output:
  # the estimated q_i,j matrix (k by k)
  return(mle_mat)
}
