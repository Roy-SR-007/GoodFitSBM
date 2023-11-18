# Get.Next.Network

# underlying model: beta SBM

# objective :: given a mixed graph G = (d,b), the routine returns new mixed graph G' in the beta-SBM fiber
# the move could be only directed, or only bidirected, or a composite of the two viz., a mixed move

#	Input::
#	d: igraph object, directed
#	b: igraph object, undirected

#	Optional Input::
#	ed.coin: vector of floats, length 3; to be used for p1.ed.recip model (we don't need it as we are only dealing with beta-SBM)
# c[1] = P(directed move); 	c[2] = P(bidirected move); c[3] = P(mixed move)
# beta.SBM.coin: a fair coin by default
# SBM.blocks: vector of integers representing the block assignment of the vertices of the undirected graph b;
# length equal to the number of vertices of b; a mandatory input for the beta.SBM model

# Output::
#  1) new.directed.graph,
#  2) new.bidirected.graph,
#  3) boolean flag trivial.move

Get.Next.Network = function(d, b, ed.coin = c(1/3, 1/3, 1/3), beta.SBM.coin = c(1/2), SBM.blocks = NULL) {

  # routine checks on the blocks
  if (is.null(SBM.blocks) || !is.vector(SBM.blocks)) {

    stop("beta.SBM model requires a non-empty vector SBM.blocks input." )
  }
  else if (length(SBM.blocks) != vcount(b)) {

    stop("Get.Next.Network error: SBM.blocks must be same length as number of vertices in b.")

  }

  # making a call to `Get_Move_beta_SBM.R`
  move = Get.Move.beta.SBM(b, blocks = SBM.blocks, coin = beta.SBM.coin)
  trivial.move = move[[3]]

  # b (undirected) minus/deleted bidirected.to.be.removed plus/added bidirected.to.be.added
  new.bidirected.graph = graph.union(graph.difference(b, move[[1]]), move[[2]])
  new.directed.graph = d

  return(list(new.directed.graph, new.bidirected.graph, trivial.move))

}
