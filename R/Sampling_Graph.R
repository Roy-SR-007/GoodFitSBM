# sample_a_move

# underlying model: beta SBM

# objective :: to sample a graph in the same fiber

# Input::
# G_current: `igraph` object which is an undirected graph with no self loop
# C: vector of block assignments of size n; block assignments varies over 1 to k

# Output::
# the graph after one random move

#' @import igraph

sample_a_move = function(C, G_current) {

  n = length(igraph::V(G_current)) # no. of vertices of the current graph

  G_comp = igraph::graph.complementer(G_current, loops=FALSE) # complement graph of the current graph

  # assign the attributes to the graph
  igraph::V(G_current)$block_asgn = C
  num_blocks = length(unique(C)) # no. of blocks

  # label edges based on the node attributes
  num_edges = length(igraph::E(G_current)) # no. of edges in the current graph
  all_edges = igraph::get.edgelist(G_current) # edge list of the current graph
  comp_edges = igraph::get.edgelist(G_comp) # edge list of the complementary graph (w.r.t to the current graph)

  # determine the type of move in the Markov basis: inter-block move (type == 1) or intra-block move (type == 2)
  # or quadratic and cubic moves switching edges along a four-cycle (type == 3)
  type = sample.int(3, size = 1)

  if(type == 1) {

    # interchanging two edges in the same block: intra-block move

    # sample a block
    s = sample.int(num_blocks, size = 1)

    # find vertices within the sampled block: vertices to add and delete for both the current graph and its complement
    to_delete = all_edges[(C[as.numeric(all_edges[ , 1])] == s) * (C[as.numeric(all_edges[ , 2])] == s) > 0, ]
    to_add = comp_edges[(C[as.numeric(comp_edges[ , 1])] == s) * (C[as.numeric(comp_edges[ , 2])] == s) > 0, ]

    # non-zero number of vertices to add and delete: sampled block has at least one edge
    if((length(to_delete) > 0) * (length(to_add) > 0)) {

      # sample an edge to add from complement graph and delete from the current graph
      delete_edge = sample.int(length(to_delete) / 2, 1) # edge to be deleted
      add_edge = sample.int(length(to_add) / 2, 1) # edge to be added

      # separating two cases: vertices to be added and deleted equals 2
      # vertices to be added and deleted not equals 2, i.e.,

      # if the sampled block has only one edge in the complement graph, then that will be added,
      # otherwise one edge will be added by sampling randomly from the complement graph
      if(length(to_add) == 2) {

        G_sample = igraph::graph.union(G_current, graph(to_add, n=n, directed = FALSE))

      }
      else {

        G_sample = igraph::graph.union(G_current, graph(to_add[add_edge, ], n = n, directed = FALSE))

      }

      # if the sampled block has only one edge in the current graph, then that will be deleted,
      # otherwise it will delete one edge by sampling randomly from the current graph
      if(length(to_delete) == 2) {

        G_sample = igraph::graph.difference(G_sample, graph(to_delete, n = n, directed = FALSE))

      }
      else {

        G_sample = igraph::graph.difference(G_sample, graph(to_delete[delete_edge, ], n = n, directed = FALSE))

      }

    }
    else {

      G_sample = G_current

    }

  }
  else if(type == 2) {

    # replace one inter edge with another between the same block pairs: inter-block move

    # sample a pair of two different blocks
    two_blocks = sample.int(num_blocks, 2)
    s = two_blocks[1]
    t = two_blocks[2]

    # check feasibility; inter edges of the current and its complementary graph
    # find edges between two fixed blocks for both the current graph and its complement graph
    inter = all_edges[((C[as.numeric(all_edges[ , 1])] == s) * (C[as.numeric(all_edges[ , 2])] == t)) + ((C[as.numeric(all_edges[ , 1])] == t) * (C[as.numeric(all_edges[ , 2])] == s)) > 0, ]
    comp_inter = comp_edges[((C[as.numeric(comp_edges[ , 1])] == s) * (C[as.numeric(comp_edges[ , 2])] == t)) + ((C[as.numeric(comp_edges[ , 1])] == t) * (C[as.numeric(comp_edges[ , 2])] == s)) > 0, ]

    # non zero inter edges: sample block has at least one edge
    if((length(inter) > 0) * (length(comp_inter) > 0)) {

      # inter edges to add and delete
      # sample an edge to add from complement graph and delete from the current graph
      delete_edge = sample.int(length(inter) / 2, 1)
      add_edge = sample.int(length(comp_inter) / 2, 1)

      # separating two cases: inter edges of the current graph and its complement graph equals 2
      # inter edges of the same not equals 2, i.e.,

      # if the sampled blocks have only one between edge in the complement graph, then that will be added,
      # otherwise it will add one between edge by sampling randomly from the complement graph
      if(length(comp_inter) == 2) {

        G_sample = igraph::graph.union(G_current, graph(comp_inter, n = n, directed = FALSE))

      }
      else {

        G_sample = igraph::graph.union(G_current, graph(comp_inter[add_edge, ],n = n, directed = FALSE))

      }

      # if the sampled blocks have only one between edge in the current graph, then that will be added,
      # otherwise it will add one between edge by sampling randomly from the current graph
      if(length(inter) == 2) {

        G_sample = igraph::graph.difference(G_sample, graph(inter, n = n, directed = FALSE))

      }
      else{

        G_sample = igraph::graph.difference(G_sample, graph(inter[delete_edge, ], n = n, directed = FALSE))

      }

    }
    else {

      G_sample = G_current

    }

  }
  else if(type == 3) {

    # four cycles with quadratic and cubic moves; making a call to `Get_Next_Network.R` routine
    list_g = Get.Next.Network(graph.empty(), G_current, ed.coin = c(0, 1, 0), SBM.blocks = C) # the `ed.coin` argument is not required in case of beta-SBM
    G_sample = list_g[[2]] # the bidirected graph from the `Get.Next.Network()` routine

  }

  # returning a graph after a sample (Markov) move
  return(G_sample)
}
