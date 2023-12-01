# block_est

# objective:: estimates in case of latent block situation

# Input::
# A a n by n binary symmetric adjacency matrix representing a undirected graph where n is the no nodes in the graph
# K a numeric scalar representing no of blocks

# Output::
# cluster: a vector of size n representing block assignment for each node; values are 1 to K i.e, no of cluster

#' @import irlba
#' @importFrom stats kmeans

block_est = function (A, K)
{

  SVD = irlba::irlba(A, nu = K, nv = K)

  km = stats::kmeans(SVD$v[, 1:K], centers = K, nstart = 30, iter.max = 100)
  return(cluster = km$cluster)

}
