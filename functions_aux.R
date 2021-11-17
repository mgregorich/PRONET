########################################################
# Author: MG
# Date: 18.10.2021
# Info: Helper functions for the simulation study
######################################################

# ------------------- GENERAL -----------------------------

VecToSymMatrix <- function(diag.entry, side.entries, mat.size, byrow=T){
  # Generate symmetric matrix from vectors holding diagonal values and side entries
  # diag.entry=0; side.entries=grow
  diags = rep(diag.entry, mat.size)
  mat = diag(diags)
  mat[lower.tri(mat, diag=FALSE)] <- side.entries
  mat <- t(mat)
  mat[lower.tri(mat, diag=FALSE)] <- side.entries   
  return(as.matrix(mat))
}



# ------------------- MODELLING -----------------------------

get_error = function(yhat, y) {
  mean.hat <- apply(yhat,1, function(x) mean(x, na.rm = T))
  bias = mean.hat - y
  relbias = ((mean.hat - y)/y)*100
  var = apply(yhat,1, function(x) var(x,na.rm=T))
  
  rmse = apply(((yhat - y) ^ 2), 1, function(x) sqrt(mean(x, na.rm = T)))
  
  out <- cbind(bias, relbias, rmse, var, "mean"=mean.hat)
  return(out)
}


# ------------------- NETWORK -----------------------------

calcGraphFeatures <- function(adj, weighted=F){
  # Calculate weighted/unweighted graph-theoretical features
  if(weighted){
    cc.w <- clustcoeff(adj, weighted = T)$CC # Clustering coefficient
    out <- c(cc.w)
  }else{
    cc.uw <- clustcoeff(adj, weighted = F)$CC  
    out <- c(cc.uw)
  }

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



evalSSN <- function(eweights, msize, tseq){
  # Perform weight-based thresholding and network feature computation (wrapper function)
  adj <- VecToSymMatrix(diag.entry = 0, side.entries = unlist(eweights), mat.size = msize)
  
  graph.vars <- sapply(tseq, function(x) calcGraphFeatures(weightThresholding(adj, x)$adj))
  return(graph.vars)
}

weightThresholding <- function(adj, w=0.5){
  # Apply weight-based thresholding to adjacency matrix
  adj[adj < w] <- 0

  return(list(adj=adj, thresh=w))
} 

