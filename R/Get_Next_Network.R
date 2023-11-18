Get.Next.Network <- function(d, b, model="p1.recip.ed", ed.coin=c(1/3,1/3,1/3), nzconst.coin=c(ecount(b)/(ecount(d)+ecount(b)), ecount(d)/(ecount(d)+ecount(b))), beta.SBM.coin=c(1/2), SBM.blocks=NULL) {

  if (is.null(SBM.blocks) || !is.vector(SBM.blocks))
    stop("beta.SBM model requires a non-empty vector SBM.blocks input." )
  else if (length(SBM.blocks)!=vcount(b))
    stop("Get.Next.Network error: SBM.blocks must be same length as number of vertices in b.")

  move = Get.Move.beta.SBM(b, blocks=SBM.blocks, coin = beta.SBM.coin)
  trivial.move = move[[3]]
  #b minus bidirected.to.be.removed plus bidirected.to.be.added
  new.bidirected.graph = graph.union(graph.difference(b, move[[1]]), move[[2]])
  new.directed.graph = d
    #    print(paste("New Network"))                        #for testing
    #    print(get.edgelist(new.bidirected.graph))          #for testing

  return(list(new.directed.graph,new.bidirected.graph,trivial.move))

}
