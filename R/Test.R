
n1 = 50
n2 = 50
n3 = 50

n = n1 + n2 + n3

class = rep(c(1, 2, 3), c(n1, n2, n3))

cmat = matrix(c(30, 0.05, 0.05,
                0.05, 30, 0.05,
                0.05, 0.05, 30),
              ncol = 3,
              byrow = TRUE)
pmat = cmat / n

adj = matrix(0, n, n)

for(i in 2:n) {

  for(j in 1:(i-1)) {

    p = pmat[class[i], class[j]]
    adj[i, j] = rbinom(1, 1, p)

  }

}

adjsymm = adj + t(adj)

goftest(adjsymm, C = class, numGraphs = 100)
