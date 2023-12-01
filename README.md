
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GoodFitSBM: Monte Carlo goodness-of-fit tests for Stochastic Blockmodels (SBMs)

<!-- <img src="man/figures/logo.png" align="right" width="150"/> -->
<!-- [![CRAN_Status_Badge](https://img.shields.io/cran/v/MatchIt?color=952100)](https://cran.r-project.org/package=lmw) [![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/MatchIt?color=952100)](https://cran.r-project.org/package=lmw) -->

### NEWS

GoodFitSBM (Version 0.0.1): `GoodFitSBM` comprises functionality that
performs the *goodness-of-fit* test for an **beta-SBM** (one of the
three variants of SBMs discussed in [Karwa et
al. (2023)](https://doi.org/10.1093/jrsssb/qkad084) used for modelling
network data.

### Note

### Overview

*Stochastic blockmodels* (SBMs) contributed to the theoretical and
algorithmic developments for analyzing *network* data, which has been
facilitated by the availability of data in diverse fields like *social
sciences*, *web recommender systems*, *protein networks*, *genomics*,
and *neuroscience*. SBMs being the generalization of the Erdős-Rényi
model given by [Erdős and Rényi
(1960)](https://doi.org/10.1515/9781400841356.38) was proposed
originally in context of *social sciences* for directed and undirected
graphs, whereas now it has been vastly extended and utilized in *latent
blocks in undirected graphs*, *latent space models*, *variable degree
distribution*, *dynamically evolving networks*, etc., becoming one of
the more popular approaches to model network data in computer science,
statistics and machine learning.

[Karwa et al. (2023)](https://doi.org/10.1093/jrsssb/qkad084) addresses
a very important aspect of *model fitting* by constructing (exact)
goodness-of-fit tests (under finite-sample setting) for three variants
of SBMs respectively used for modeling network data (where model
adequacy procedures are somewhat elusive in general), viz., *Erdős-Rényi
SBM* (ER-SBM), *Additive SBM*, and *beta-SBM*, where the main idea
revolves around a *frequentist conditional goodness-of-fit test*
conditioned on a sufficient statistic ([Karwa et
al. (2023)](https://doi.org/10.1093/jrsssb/qkad084) also illustrates its
Bayesian counterpart).

### Intended Use of the package `GoodFitSBM`

Out of the three variants of the SBMs in ([Karwa et
al. (2023)](https://doi.org/10.1093/jrsssb/qkad084)) to model network
data, `GoodfitSBM` addresses goodness-of-fit test under the framework of
a **beta-SBM**. With a focus on *simple undirected*, and *unweighted*
networks (graphs) having *no self-loop*, the package comprises of four
functions viz., `get_mle()`, `goftest()`, `graphchi()`, and
`sample_a_move()` - among which `goftest()` performs the major
functionality of performing the goodness-of-fit test, in process,
obtains the value of the chi-square test statistic and the corresponding
$p-$value (using a Monte Carlo approach), after proper sampling of the
graph (a Markov move or basis) - done by the function `sample_a_move()`
under the beta-SBM framework. Computation of the chi-square test
statistic value under the beta-SBM framework is done by `graphchi()`,
which in turn requires the estimates (maximum likelihood estimation
(MLE) in our case) of the edge probabilities $q_{ij}$ between block
$B_i$ and block $B_j$ as, $\widehat{q}_{ij}$, which is done by the
function `get_mle()`.

The corresponding $p-$value obtained from `goftest()` are further
analysed to examine the extent of fit of the beta-SBM to the given
network (graph) (usual rule applies to reject the null of a good fit
when, $p \leq \alpha$ at level $0< \alpha < 1$). There are other
functions collated, each performing the required functionality in the
process, see [GoodFitSBM Github
Repo](https://github.com/Roy-SR-007/GoodFitSBM/tree/master/R) for
details.

### The beta-SBM and Related Goodness-of-Fit (GoF) Framework

In this section, we outline the theoretical part on which the entire
package `GoodFitSBM` is based.

#### Stochastic Blockmodels (SBMs)

Consider $G$ to be a graph on $n$ nodes, $g$ being its realization. It
is assumed that, all graphs are *unweighted* and *undirected*, and there
are *no self-loops*. Therefore, the graph $G$ can also be referred to by
its $n\times n$ adjacency matrix, $\mathbf{g}$, where $g_{uv} = 1$, if
there is an edge between node $u$ and $v$ in the graph $\mathbf{g}$, and
$0$ otherwise.

An *exponential family random graph model (ERGM)* assumes that the
probability of observing a graph $\mathbf{g}$ depends only on a vector
of sufficient statistics, $T(\mathbf{g})$, i.e., the probability of
observing a given network $G=\mathbf{g}$ is,

$$
\mathbb{P}_{\theta}(G=\mathbf{g}) = \frac{exp{\langle T(\mathbf{g}), \theta\rangle}}{\psi(\theta)}
$$ where
$\psi(\theta) = \sum_{\mathbf{g}}exp(\langle T(\mathbf{g}), \theta\rangle)$
is the normalizing constant, $\theta \in \Theta$ is a vector of natural
parameters, and $T(\mathbf{g})$ is the vector of minimal sufficient
statistics for the underlying model. Let $p_{uv}$ be the probability of
an edge between node $u$ and node $v$, where it is assumed that,
$g_{uv}\sim \mbox{Bernoulli}(p_{uv})$.

A *stochastic blockmodel (SBM)* postulates that the nodes are
partitioned into $k$ blocks and the probability $p_{uv}$ depends on
their block membership, where ([Karwa et
al. (2023)](https://doi.org/10.1093/jrsssb/qkad084)) considers three
different *log-linear* parametrizations of $p_{uv}$, one of which yields
the beta-SBM.

#### The beta-SBM

Beta-SBM is exponential family version of the degree-corrected
stochastic blockmodels suggested by [Karrer and Newman
(2011)](https://link.aps.org/doi/10.1103/PhysRevE.83.016107). The
log-linear model for the edge probabilities is parametrized by
$\binom{k+1}{2}$ block parameters $\alpha_{z(u)z(v)}$ and $n$
node-specific parameters $\beta_{u}$ for $1\leq u \leq n$, with the
log-odds of the probability of an edge $\{u,v\}$ given by,

$$
\log\bigg(\frac{p_{uv}}{1 - p_{uv}}\bigg) = \alpha_{z(u)z(v)} + \beta_{u} + \beta_{v}
$$ When $\mathcal{Z}$ (the block assignments) are known, the beta-SBM is
an exponential family model with natural parameter vector
$\theta = (\beta, \alpha)$, where $\beta = (\beta_1, \ldots, \beta_{n})$
and $\alpha$ is the vector of the upper diagonal elements of the
$k\times k$ matrix $((\alpha_{i,j}))$. The natural parameter space is
$\Theta \equiv \mathbb{R}^n \times \mathbb{R}^{\binom{k+1}{2}}$. The
vector of sufficient statistics $T_{\beta}$ contains the degree of each
node $i\in [n]$ and the number of edges between each pair of blocks
$B_{i}$ and $B_{j}$ for $1\leq i < j \leq k$.

#### Goodness-of-Fit Test for the beta-SBM

Let us define the goodness-of-fit test as,

$$
H_0: G \sim \mathbb{P}_{\theta}(G|Z)
$$ against general alternatives, where $\mathbb{P}_{\theta}(G|Z)$ is a
variant of the SBM with a generic parameter vector $\theta$ and block
assignment $Z$ (it is assumed that the number of blocks $k$ is fixed and
known).

Considering $T(\mathbf{g})$ to be the vector of sufficient statistics in
$\mathbb{P}_{\theta}(G|Z)$, define a subset of the sample space as,
$F_{u} := \{\mathbf{g}:T(\mathbf{g}) = u\}$, which is the known as the
*fiber* of $u$ under the given exponential family model. Now the usual
chi-square statistic (for the goodness-of-fit test) is not constant on
the above defined fibers of the beta-SBM variant. Thus, analogous to the
beta-model in [Petrovic et al. (2010)](https://arxiv.org/abs/0909.0073),
we use the Pearson’s chi-square statistic,

$$
GoF_{Z}(\mathbf{g}) = \chi^{2}_{\mbox{beta-SBM}}(\mathbf{g}, Z) = \sum_{1\leq u < v \leq n}\frac{(g_{uv} - \widehat{g}_{uv})^2}{\widehat{g}_{uv}}
$$ where,
$\widehat{g}_{uv} = exp(\widehat{\alpha}_{z(u)z(v)} + \widehat{\beta}_{u} + \widehat{\beta}_{v})/( exp(\widehat{\alpha}_{z(u)z(v)} + \widehat{\beta}_{u} + \widehat{\beta}_{v}))$,
where the MLEs $\widehat{\alpha}$ and $\widehat{\beta}$ are associated
with the MLEs of the edge probabilities $q_{uv}$, as obtained in the
`get_mle()` function. Large values of $\chi^{2}_{\mbox{beta-SBM}}$
correspond to departure from the null $(H_0)$ as stated above.

- For the algorithm of sampling of the graph (as done by
  `sample_a_move()`) - or a Markov move (basis) - and other details, see
  Section 5.3 of ([Karwa et
  al. (2023)](https://doi.org/10.1093/jrsssb/qkad084)).

### Use and Related Functionalities

In this section, we highlight the installation and different
functionalities included in `GoodFitSBM` along with its implementation
on some test cases and a real-life dataset.

#### Installation

To install and load-up the (development version 0.0.1) package
`GoodFitSBM` from [GoodFitSBM Github
Repo](https://github.com/Roy-SR-007/GoodFitSBM), run the following
commands.

``` r1
# install.packages("devtools")
# install.packages("remotes")
remotes::install_github("Roy-SR-007/GoodFitSBM")

library(GoodFitSBM)
```

#### Example 1: Sampling a Graph using `sample_a_move()`

Here we consider sampling (Markov move) a graph (under the beta-SBM
framework) with a total of $n=15$ nodes, with $k = 3$ blocks of sizes
$5$ each.

``` r2
library(igraph)
RNGkind(sample.kind = "Rounding")
set.seed(1729)

# We model a network with 3 even classes
n1 = 5
n2 = 5
n3 = 5

# Generating block assignments for each of the nodes
n = n1 + n2 + n3
class = rep(c(1, 2, 3), c(n1, n2, n3))

# Generating the adjacency matrix of the network
# Generate the matrix of connection probabilities
cmat = matrix(
  c(
    30, 0.05, 0.05,
    0.05, 30, 0.05,
    0.05, 0.05, 30
  ),
  ncol = 3,
  byrow = TRUE
)
pmat = cmat / n

# Creating the n x n adjacency matrix
adj <- matrix(0, n, n)
for (i in 2:n) {
  for (j in 1:(i - 1)) {
    p = pmat[class[i], class[j]] # We find the probability of connection with the weights
    adj[i, j] = rbinom(1, 1, p) # We include the edge with probability p
  }
}

adjsymm = adj + t(adj)

# graph from the adjacency matrix
G = igraph::graph_from_adjacency_matrix(adjsymm, mode = "undirected", weighted = NULL)

# plotting the current graph
plot(G, main = "The current graph")

# sampling a Markov move for the beta SBM
G_sample = sample_a_move(class, G)

# plotting the sampled graph
plot(G_sample, main = "The sampled graph after one Markov move for beta SBM")
```
