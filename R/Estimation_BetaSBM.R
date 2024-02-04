#'
#' @title Maximum Likelihood Estimation of edge probabilities between blocks of a graph
#'
#' @description `get_mle_BetaSBM` obtains MLE for the probability of edges between blocks in a graph, used in calculating the goodness-of-fit test statistic for the beta-SBM (Karwa et al. (2023))
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
#' @seealso [goftest_BetaSBM()] performs the goodness-of-fit test for the beta-SBM, where the MLE of the edge probabilities are required
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
#' get_mle_BetaSBM(G, class)
#'
#' @references
#' Karwa et al. (2023). "Monte Carlo goodness-of-fit tests for degree corrected and related stochastic blockmodels",
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology},
#' <https://doi.org/10.1093/jrsssb/qkad084>

get_mle_BetaSBM = function(G, C) {

  # get_mle_BetaSBM

  # underlying model: beta-SBM

  # objective :: calculating the MLE of the probability of edges between blocks

  # Input::
  # G: `igraph` object which is an undirected graph with no self loop
  # C: numeric vector of size n of block assignment, from 1 to k

  # Output::
  # A MLE matrix with the MLE of the probability of edges between pairs of blocks

  # collapse to k*k
  k = length(unique(C)) # no. of blocks
  n = length(C) # no. of nodes
  edges = igraph::get.edgelist(G) # edges of the graph G
  A = igraph::get.adjacency(G , type = "both") # adjacency matrix


  table_slice = array(0, c(n, n, k))
  start_table = array(0, c(n, n, k))

  for(idyad in 1:k) {

    table_slice[ , , idyad] = as.matrix(A * ((C == idyad) %*% t(rep(1, n))))
    start_table[ , , idyad] = ((C == idyad) %*% t(rep(1, n)))

  }

  fm = stats::loglin(table_slice, list(c(3)), fit=TRUE, start = start_table) # log-linear model for the cell-probabilities
  largemle = fm$fit

  mleMatr = matrix(0, nrow = n, ncol = n)
  for(i in 1:n) {

    for(j in 1:n) {

      mleMatr[i, j] = max(sum(largemle[i, j, ]), sum(largemle[j, i, ]))

    }

  }

  return(mleMatr)

}
