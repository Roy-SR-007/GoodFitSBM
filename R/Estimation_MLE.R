# get_mle

# underlying model: beta SBM

# objective :: calculating the MLE of the probability of edges between blocks

# Input::
# G: `igraph` object which is an undirected graph with no self loop
# C: numeric vector of size n of block assignment, from 1 to k

# Output::
# A MLE matrix with the MLE of the probability of edges between pairs of blocks

get_mle = function(G, C) {

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

  fm = loglin(table_slice, list(c(3)), fit=TRUE, start = start_table) # log-linear model for the cell-probabilities
  largemle = fm$fit

  mleMatr = matrix(0, nrow = n, ncol = n)
  for(i in 1:n) {

    for(j in 1:n) {

      mleMatr[i, j] = max(sum(largemle[i, j, ]), sum(largemle[j, i, ]))

    }

  }

  return(mleMatr)

}
