########################################################
# Author: MG
# Date: 18.10.2021
# Info: Helper functions for the simulation study
#########################################################


calcGraphMetrics <- function(adj){
  # Calculate graph features from correlation matrix
  cc.uw <- clustcoeff(adj, weighted = F)$CC  # Clustering coefficient
  cc.w <- clustcoeff(adj, weighted = T)$CC
  cpl.uw <- pathlengths(adj, weighted = F)
  cpl.w <- pathlengths(adj, weighted = T)
  
  net <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = F)
  cl <- clusters(net)
  mod <- modularity(net, membership = cl$membership, weights=E(net)$weight)
  ass <- assortativity.degree(net, directed=F)
  ev <- eigen_centrality(net, weights = E(net)$weight)
  eigen.score <- ev$value
  
  out <- c(cc.uw, cc.w, cpl.uw$ASPL, cpl.w$ASPL, cpl.uw$diameter, cpl.w$diameter, mod=mod, ass=ass, eigen=eigen.score)
  colnames(out) <- c("cc.uw","cc.w","cpl.uw","cpl.w","dia.uw","dia.w", "mod", "ass", "eigen") 
  return(out)
}


densityThresholding <- function(adj, d=0.5){
  # Apply density-based thresholding to adjacency matrix
  V= ncol(adj)
  E.vals=sort(adj[upper.tri(adj)],decreasing = T)
  E.poss <- (V*(V-1))/2
  E.d <- round(E.poss*d,0)
  min.ew <- min(E.vals[1:E.d])
  adj[adj< min.ew] <- 0
  
  return(list(adj=as.matrix(adj), thresh=d))
} 

VecToSymMatrix <- function(diag.entry, side.entries, byrow=T){
  # Generate symmetric matrix from vectors holding diagonal values and side entries
  mat=diag(diag.entry)
  mat[lower.tri(mat, diag=FALSE)] <- -ox
  mat <- t(mat)
  mat[lower.tri(mat, diag=FALSE)] <- -ox   
  return(mat)
}