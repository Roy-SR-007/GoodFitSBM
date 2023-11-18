Get.Directed.Piece <- function(d){
  # d = directed part of G.
  # pick a random subset E of edges of d and randomly shuffle it
  # (i.e., E = random sample from d of random size):
  if (ecount(d)==2) {# avoid unwanted behaviour of sample function
    random.subset.of.d=get.edges(d,1:2)
    subset.size=2
  }
  else if (ecount(d)>2){
    subset.size = sample(2:ecount(d) ,1) #this is a random integer
    random.edge.indices = sample(1:(ecount(d)),subset.size)
    random.subset.of.d = get.edges(d,random.edge.indices)
  }
  else return(NULL)
  # randomly partition E,
  # and for every part E_i, call Bipartite.Walk(E_i)
  # and merge the edges.to.add_i from each of the partitions into a big set edges.to.add
  number.of.partitions = sample(1:(floor(subset.size/2)), 1)
  # initialize where to store the pieces of the walk:
  edges.to.add = c()
  edges.to.remove = c()
  more.edges = c()
  num.edges.left =subset.size
  s=1 #index
  while(num.edges.left>1) {
    if (num.edges.left==2) k=2 #avoid unwanted behaviour of sample function
    else k = sample(2:num.edges.left,1) #size of current part.
    if (num.edges.left-k == 1) k=k+1 #E's assumption on not leaving out that last edge hanging.
    more.edges=Bipartite.Walk(random.subset.of.d[s:(s+k-1),])
    if (is.null(more.edges)) return(NULL)
    else edges.to.add = c(edges.to.add,more.edges )
    num.edges.left=num.edges.left-k
    s=s+k
  }
  # edges.to.remove has to be in the same format as edges.to.add, so do this:
  if ( !is.null(edges.to.add)) as.vector(t(random.subset.of.d)) -> edges.to.remove
  return(list(edges.to.remove,edges.to.add))
}
#######################################################################
#######################################################################
Get.Bidirected.Piece <- function(b) {
  ## THIS function computes bidirected move ONLY without checks for conflicts.
  # this calls Bipartite.Walk but first checks if edges are a matching?
  # Randomly direct the entire bidirected graph and call Get.Directed.Piece

  if (ecount(b) < 2)
    return(NULL)
  b.directed = as.arbitrary.directed(b)
  return(Get.Directed.Piece(b.directed))
}

as.arbitrary.directed <- function(b) {
  # Create a directed graph out of the edges of b
  b.decr = graph(t(get.edges(b, 1:ecount(b))))
  # Pick a random integer from 0 to #edges in b
  num.edges.to.reverse = sample(0:ecount(b), 1)
  # Direct the first num.edges.to.reverse edges in one way and the others the other way
  if (num.edges.to.reverse==0) {
    b.directed = b.decr
  }else{
    random.edge.indices = sample(1:ecount(b), num.edges.to.reverse)
    b.subset.decr = graph(t(get.edges(b, random.edge.indices))) #get.edges and get.edgelist direct edges in a different order somehow!
    el = get.edgelist(b.subset.decr, names = FALSE) ## magically swap cols to reverse direction
    b.subset.incr = graph(rbind(el[, 2], el[, 1]))
    # make the directed graph out of:
    # (reversed.edges.of.b being directed in the decreasing order) union (remaining edges in incr.order):
    b.directed = graph.union(graph.difference(b.decr, b.subset.decr), b.subset.incr)
  }
  return(b.directed)
}
