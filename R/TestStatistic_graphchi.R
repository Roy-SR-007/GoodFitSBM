# graphchi

# underlying model: beta SBM

# objective :: calculate the chi-square statistics given the graph and the MLE

# Input::
# G: igraph` object which is an undirected graph with no self loop
# C: numeric vector of size n of block assignment, from 1 to k
# p_mle: matrix of MLE table

# Output::
# the value for chi-sq statistics

graphchi = function(G, C, p_mle) {

  # getting graph information
  n = length(V(G))
  k = length(unique(C))
  C = as.vector(C)
  A = igraph::get.adjacency(G, type = "both")

  blockchi = matrix(0, nrow = k, ncol = k)

  # observed no. of neighbors of an element, in each block
  degseq = matrix(0, nrow = n, ncol = k)

  # expected no. of neighbors of an element, in each block
  exp_mat = matrix(0, nrow = n, ncol = k)

  for(inode in 1 : n) {

    row = as.numeric(A[inode, ])

    for(iblock in 1 : k) {

      degseq[inode, iblock] = sum(row[C == iblock])
      ni = sum(C == iblock) # no. of nodes in the ith block
      exp_mat[inode, iblock] = ni * mean(p_mle[C == iblock, inode])

    }

  }

  # returning the value of the chi-square test statistics
  return(sum((degseq[exp_mat != 0] - exp_mat[exp_mat != 0]) ^ 2 / exp_mat[exp_mat != 0]))

}
