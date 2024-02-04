#'
#' @title Computation of the chi-square test statistic for goodness-of-fit under ERSBM
#' @description `graphchi_ERSBM` obtains the value of the chi-square test statistic required for the goodness-of-fit of a ERSBM (Karwa et al. (2023))
#'
#' @param G an igraph object which is an undirected graph with no self loop
#' @param C a positive integer vector of size n for block assignments of each node; from 1 to K (no of blocks)
#' @param p_mle a matrix with the MLE estimates of the edge probabilities
#'
#' @return A numeric value
#' \item{teststat_val}{The value of the chi-square test statistic}
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
#'
#' @export
#'
#' @seealso [goftest_ERSBM()] performs the goodness-of-fit test for the ERSBM, where the values of the chi-square test statistics are required
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
#' p.hat = get_mle_ERSBM(G, class)
#'
#' # chi-square test statistic values
#' graphchi_ERSBM(G, class, p.hat)
graphchi_ERSBM <- function(G, C, p_mle) {
  # Input: G: igraph object which is an undirected graph and has no self loop
  #        C: numeric vector of size n of block assignment; from 1 to k
  #        p_mle: k*k matrix of MLE table

  # Getting graph information
  n <- length(C)
  k <- length(unique(C))
  A <- igraph::get.adjacency(G, type = "both")

  # Observed no of neighbors of a element in each block
  degseq <- matrix(0, nrow = n, ncol = k)

  # Expected no of neighbors of a element in each block
  exp_mat <- matrix(0, nrow = n, ncol = k)
  for (inode in 1:n) {
    row <- as.numeric(A[inode, ])
    for (iblock in 1:k) {
      degseq[inode, iblock] <- sum(row[C == iblock])
      ni <- sum(C == iblock) # no of nodes in ith block
      exp_mat[inode, iblock] <- ni * p_mle[iblock, C[inode]]
    }
  }

  # Output:
  # the value for chi-sq statistics
  return(sum((degseq[exp_mat != 0] - exp_mat[exp_mat != 0])^2 / exp_mat[exp_mat != 0]))
}
