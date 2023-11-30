# Usual Testing of goodness-of-fit test in case of beta SBM, for modeling network data

#' @import igraph

# Test 1

set.seed(1729)

# We model a network with 3 even classes
n1 = 50
n2 = 50
n3 = 50

# Generating block assignments for each of the nodes
n = n1 + n2 + n3
class = rep(c(1, 2, 3), c(n1, n2, n3))

# Generating the adjacency matrix of the network
# Generate the matrix of connection probabilities
cmat = matrix(c(30, 0.05, 0.05,
                0.05, 30, 0.05,
                0.05, 0.05, 30),
              ncol = 3,
              byrow = TRUE)

pmat = cmat / n

# Creating the n x n adjacency matrix
adj = matrix(0, n, n)

for(i in 2:n) {

  for(j in 1:(i-1)) {

    p = pmat[class[i], class[j]]
    adj[i, j] = rbinom(1, 1, p)

  }

}

adjsymm = adj + t(adj)

# When class assignment is known
out = goftest(adjsymm, C = class, numGraphs = 100)

chi_sq_seq = out$statistic
pvalue = out$p.value
print(pvalue)

# Plotting histogram of the sequence of the test statistics
hist(chi_sq_seq, 20, xlab = "chi-square test statistics", main = NULL)
abline(v = chi_sq_seq[1], col = "red", lwd = 5) # adding test statistic on the observed network


# Test 2

set.seed(1729)

# We model a network with 3 classes
n1 = 30
n2 = 20
n3 = 50

# Generating block assignments for each of the nodes
n = n1 + n2 + n3
class = rep(c(1, 2, 3), c(n1, n2, n3))

# Generating the adjacency matrix of the network
# Generate the matrix of connection probabilities
cmat = matrix(c(30, 0.05, 0.05,
                0.05, 30, 0.05,
                0.05, 0.05, 30),
              ncol = 3,
              byrow = TRUE)

pmat = cmat / n

# Creating the n x n adjacency matrix
adj = matrix(0, n, n)

for(i in 2:n) {

  for(j in 1:(i-1)) {

    p = pmat[class[i], class[j]]
    adj[i, j] = rbinom(1, 1, p)

  }

}

adjsymm = adj + t(adj)

# When class assignment is known
out = goftest(adjsymm, C = class, numGraphs = 100)

chi_sq_seq = out$statistic
pvalue = out$p.value
print(pvalue)

# Plotting histogram of the sequence of the test statistics
hist(chi_sq_seq, 20, xlab = "chi-square test statistics", main = NULL)
abline(v = chi_sq_seq[1], col = "red", lwd = 5) # adding test statistic on the observed network
