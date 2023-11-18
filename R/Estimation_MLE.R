get_mle = function(G,C,modelname){
  # collapse to k*k
  k = length(unique(C));
  n = length(C);
  edges=get.edgelist(G);
  A=get.adjacency(G,type="both");

  #A=as_adjacency_matrix(G,type='both');
  table_slice = array(0, c(n, n, k));
  start_table = array(0, c(n, n, k));
  for(idyad in 1:k){
    table_slice[,,idyad]=as.matrix(A*((C==idyad)%*%t(rep(1,n))));
    start_table[,,idyad]=((C==idyad)%*%t(rep(1,n)));
  }
  fm <- loglin(table_slice, list(c(3)), fit=TRUE, start=start_table);
  largemle = fm$fit;
  mleMatr = matrix(0,nrow=n,ncol=n);
  for(i in 1:n){
    for(j in 1:n){
      mleMatr[i,j]=max(sum(largemle[i,j,]),sum(largemle[j,i,]));
    }
  }
}
