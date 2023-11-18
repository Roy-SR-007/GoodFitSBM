graphchi = function(G,C,p_mle,modelname='betablock'){
  # Calculate the chi-square statistics given the graph and the MLE
  # Input: G: igraph object of a graph
  #        C: block assignment
  #        p_mle: k*k matrix of MLE table
  # Output: the value for chi-sq statistics

  n=length(V(G)); k=length(unique(C));C=as.vector(C);
  A=get.adjacency(G,type="both");

  #A=as_adjacency_matrix(G,type='both')
  blockchi=matrix(0,nrow=k,ncol=k);
  degseq = matrix(0,nrow=n,ncol=k);
  exp_mat = matrix(0,nrow=n,ncol=k);
  for(inode in 1:n){
    row = as.numeric(A[inode,]);
    for(iblock in 1:k){
      degseq[inode,iblock]=sum(row[C==iblock]);
      ni=sum(C==iblock);
      exp_mat[inode,iblock] = ni*mean(p_mle[C==iblock,inode]);
    }
  }
  return(sum((degseq[exp_mat!=0]-exp_mat[exp_mat!=0])^2/exp_mat[exp_mat!=0]));
  #return(sum(blockchi))
}
