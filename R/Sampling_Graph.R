# Exact testing for block models

sample_a_move = function(C,G_current,modelname='betablock') {
  # Sample a graph in the same fiber
  # Input:
  # G_current: igraph object G_obs
  # C: vector of class assignment
  # modelname: 'betablock' or 'sbm'
  # Output:
  # the graph after one random move
  n=length(V(G_current));
  G_comp = graph.complementer(G_current, loops=FALSE)
  if(modelname=='betablock'){
    # Assign the attributes to the graph
    V(G_current)$block_asgn=C
    num_blocks = length(unique(C))
    # Label edges based on the node attributes
    num_edges = length(E(G_current));
    all_edges = get.edgelist(G_current);
    comp_edges = get.edgelist(G_comp);
    # Determine the type
    type = sample.int(3,size=1);
    if(type==1){
      # interchanging two edges in the same block
      # sample a block
      s = sample.int(num_blocks,1);
      # find vertexs within the block
      to_delete = all_edges[(C[as.numeric(all_edges[,1])]==s)*(C[as.numeric(all_edges[,2])]==s)>0,];
      to_add = comp_edges[(C[as.numeric(comp_edges[,1])]==s)*(C[as.numeric(comp_edges[,2])]==s)>0,];
      if((length(to_delete)>0)*(length(to_add)>0)){
        delete_edge = sample.int(length(to_delete)/2,1);
        add_edge = sample.int(length(to_add)/2,1);
        if(length(to_add)==2){
          G_sample = graph.union(G_current,graph(to_add,n=n,directed=FALSE))}
        else{
          G_sample = graph.union(G_current,graph(to_add[add_edge,],n=n,directed=FALSE))
        }
        if(length(to_delete)==2){
          G_sample=graph.difference(G_sample,graph(to_delete,n=n,directed=FALSE))
        }
        else{
          G_sample=graph.difference(G_sample,graph(to_delete[delete_edge,],n=n,directed=FALSE))
        }
      }
      else{
        G_sample=G_current;
      }
    }else if(type==2){
      # replace one inter edge with another between the same block pairs
      two_blocks = sample.int(num_blocks,2);
      s = two_blocks[1];t = two_blocks[2];
      # Check feasibility
      inter = all_edges[((C[as.numeric(all_edges[,1])]==s)*(C[as.numeric(all_edges[,2])]==t))+((C[as.numeric(all_edges[,1])]==t)*(C[as.numeric(all_edges[,2])]==s))>0,];
      comp_inter = comp_edges[((C[as.numeric(comp_edges[,1])]==s)*(C[as.numeric(comp_edges[,2])]==t))+((C[as.numeric(comp_edges[,1])]==t)*(C[as.numeric(comp_edges[,2])]==s))>0,];
      if((length(inter)>0)*(length(comp_inter)>0)){
        delete_edge = sample.int(length(inter)/2,1);
        add_edge = sample.int(length(comp_inter)/2,1);
        if(length(comp_inter)==2){
          G_sample = graph.union(G_current,graph(comp_inter,n=n,directed=FALSE))
        }
        else{
          G_sample = graph.union(G_current,graph(comp_inter[add_edge,],n=n,directed=FALSE))
        }
        if(length(inter)==2){
          G_sample=graph.difference(G_sample,graph(inter,n=n,directed=FALSE))
        }else{
          G_sample=graph.difference(G_sample,graph(inter[delete_edge,],n=n,directed=FALSE))
        }
      }else{
        G_sample=G_current;
      }
    }else if(type==3){
      # four cycles
      list_g = Get.Next.Network(graph.empty(),G_current, ed.coin=c(0,1,0));
      G_sample = list_g[[2]];
    }

    return(G_sample);
  }
}
