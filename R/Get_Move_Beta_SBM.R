# Get.Move.beta.SBM

# underlying model: beta SBM

# objective :: returning a random move, not necessarily primitive that is applicable to the observed
# network (graph) G that is guaranteed to move to a network in the beta-SBM fiber.

# Input::
# g: undirected graph, as an undirected igraph object
# blocks: a vector designating the block assignments of the nodes; blocks[i] = j if vertex i is assigned to block j
# coin: controls

# Output::
# a list of two undirected igraph objects representing the move (edges to remove and add) and
# boolean flag indicating whether the move is empty

#' @import igraph
#' @include Get_Induced_Subgraph.R
#' @include Get_Bidirected_Move.R
#' @include Get_Between_Blocks_Move_beta_SBM.R
#' @include Get_Within_Blocks_beta_SBM.R

Get.Move.beta.SBM = function(g, blocks, coin = c(1/2)) {

  # iterate through all the blocks generating a random bidirected move within
  # each block and a bipartite move between blocks.
  if(is.null(blocks)) {

    print("Error: blocks parameter cannot be empty in Get.Move.beta.SBM.")

  }

  n = igraph::vcount(g) # no. of vertices of the graph g
  k = max(blocks)

  move = list(igraph::graph.empty(n, directed=FALSE), igraph::graph.empty(n, directed=FALSE), TRUE)

  v.block = list()
  g.block = list()

  for (i in 1:k) {

    # The subgraphs within the blocks

    v.block[[i]] = which(blocks == i)
    g.block[[i]] = Get.Induced.Subgraph(g, v.block[[i]])

  }

  coin.value = runif(1)
  if (coin.value <= coin[1]) {

    # this part of the code allows for moves containing the switch ab, cd<->ac, bd where block[a] = block[b] = block[c] != block[d]
    # to be produced. These moves are both valid and necessary as such moves preserve everyone's degree in the full graph,
    # allowing degrees of individual vertices to change within a block and between two blocks, and do not change the number of edges
    # within block i or between block i and block j.

    small.move.coin.value = runif(1)
    if(small.move.coin.value < 0.5) {

      # perform a small move guaranteed to generate a network in a different sub-fiber with different within-block vertex degrees

      indices = sample(1:k, 2)
      i = indices[1]
      j = indices[2]

      if (length(v.block[[i]] >= 3) && length(v.block[[j]]) > 1) {

        v.included = c(sample(v.block[[i]], 3), sample(v.block[[j]], 1))

      }
      else if (length(v.block[[i]] >= 3) && length(v.block[[j]] == 1)) {

        v.included = c(sample(v.block[[i]], 3), v.block[[j]])

      }
      else {

        return (move)

      }

    }
    else{

      r = sample(2:k, 1)
      included.blocks = sample(1:k, r)
      v.included = c()
      for (i in included.blocks) {

        v.included = c(v.included, v.block[[i]])

      }

    }

    g.subgraph = Get.Induced.Subgraph(g, v.included)
    proposed.move = list(igraph::graph.empty(n), igraph::graph.empty(n, directed = FALSE))

    count = 0
    while (igraph::ecount(proposed.move[[1]]) == 0 && count < 50) {

      proposed.move = Get.Bidirected.Move(graph.empty(vcount(g.subgraph)), g.subgraph)
      count = count + 1

    }

    # error checking; making sure we are removing same total number of edges as we are adding
    if (igraph::ecount(proposed.move[[1]]) != igraph::ecount(proposed.move[[2]])) {

      return(list(igraph::graph.empty(n, directed = FALSE), igraph::graph.empty(n, directed = FALSE), TRUE))

    }

    graph.to.remove = igraph::graph.difference(proposed.move[[1]], proposed.move[[2]])
    graph.to.add = igraph::graph.difference(proposed.move[[2]], proposed.move[[1]])

    if (igraph::ecount(graph.to.remove) == 0 && igraph::ecount(graph.to.add) == 0) {

      return(list(igraph::graph.empty(n, directed = FALSE), igraph::graph.empty(n, directed = FALSE), TRUE))

    }

    for (i in 1:k) {

      # check that number of edges within block i remain constant
      if (igraph::ecount(Get.Induced.Subgraph(graph.to.remove, v.block[[i]])) != igraph::ecount(Get.Induced.Subgraph(graph.to.add, v.block[[i]]))) {

        return(list(igraph::graph.empty(n, directed = FALSE), igraph::graph.empty(n, directed = FALSE), TRUE))

      }

      if (i < k) {

        for (j in (i + 1):k) {

          # check that number of edges between block i and block j remain constant
          g.full.i.j = Get.Induced.Subgraph(g, c(v.block[[1]], v.block[[2]]))
          g.between.i.j = igraph::graph.difference(g.full.i.j, igraph::graph.union(g.block[[i]], g.block[[j]]))
          num.g.between.i.j.edges.to.remove = length(which(get.edge.ids(g.between.i.j, as.vector(t(get.edgelist(graph.to.remove)))) > 0))

          num.g.between.i.j.edges.to.add = 0
          edges.to.add = get.edgelist(graph.to.add)

          for (l in 1:igraph::ecount(graph.to.add)) {

            if ((is.element(edges.to.add[l, 1], v.block[[i]]) && is.element(edges.to.add[l, 2], v.block[[j]])) || (is.element(edges.to.add[l, 1], v.block[[j]]) && is.element(edges.to.add[l, 2], v.block[[i]]))) {

              num.g.between.i.j.edges.to.add = num.g.between.i.j.edges.to.add + 1

            }

          }

          if (num.g.between.i.j.edges.to.add != num.g.between.i.j.edges.to.remove) {

            return(list(igraph::graph.empty(n, directed = FALSE), igraph::graph.empty(n, directed = FALSE), TRUE))

          }

        }

      }

    }

    move[[1]] = proposed.move[[1]] # graph.to.remove
    move[[2]] = proposed.move[[2]] # graph.to.add

    if (igraph::ecount(move[[1]]) > 0 && igraph::ecount(igraph::graph.difference(move[[1]], move[[2]])) > 0 && igraph::ecount(igraph::graph.difference(move[[2]], move[[1]]))>0) {

      move[[3]] = FALSE
    }

  }
  else {

    for (i in 1:k) {

      # produce a move within block i
      move.within.i = Get.Within.Blocks.Move.beta.SBM(g.block[[i]])
      move[[1]] = igraph::graph.union(move.within.i[[1]], move[[1]])
      move[[2]] = igraph::graph.union(move.within.i[[2]], move[[2]])

      if (i < k) {

        for (j in (i + 1):k) {

          # produce a move between block i and j
          g.full.i.j = Get.Induced.Subgraph(g, c(v.block[[1]], v.block[[2]]))
          g.between.i.j = igraph::graph.difference(g.full.i.j, igraph::graph.union(g.block[[i]], g.block[[j]]))
          move.between.i.j = Get.Between.Blocks.Move.beta.SBM(g.between.i.j)
          move[[1]] = igraph::graph.union(move.between.i.j[[1]], move[[1]])
          move[[2]] = igraph::graph.union(move.between.i.j[[2]], move[[2]])

        }

      }

    }

    # error checking making sure we are removing same number of edges as we are adding
    if (igraph::ecount(move[[1]]) != igraph::ecount(move[[2]])) {

      return(list(igraph::graph.empty(n, directed = FALSE), igraph::graph.empty(n, directed = FALSE), TRUE))

    }

    if (igraph::ecount(move[[1]]) > 0 && igraph::ecount(igraph::graph.difference(move[[1]], move[[2]])) > 0 && igraph::ecount(igraph::graph.difference(move[[2]], move[[1]])) > 0) {

      move[[3]] = FALSE

    }

  }

  return (move)

}
