########################################################
# Author: MG
# Date: 18.10.2021
# Info: Helper functions for the simulation study
#########################################################


# ------------------- FUNCTIONS -----------------------------

calcGraphFeatures <- function(adj, weighted=T){
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

evalLM <- function(x, dt){
  # Perform univariable linear regression and extract coeffs, rmse and adjusted r2
  fit.lm <- lm(as.formula(paste0("Y ~ ",x)), data=dt)
  sum.lm <- summary(fit.lm)
  rmse <-  sqrt(mean((data.gvars$Y-fit.lm$fitted.values)^2))
  
  out <- list("intercept"=as.numeric(fit.lm$coefficients[1]), "slope"=as.numeric(fit.lm$coefficients[2]),
              "adjr2"=sum.lm$adj.r.squared, "rmse"=rmse)
  return(out)
}

evalSSN <- function(eweights, msize){
  # Perform weight-based thresholding and network feature computation (wrapper function)
  adj <- VecToSymMatrix(diag.entry = 0, side.entries = unlist(eweights), mat.size = msize)
  
  tseq <- seq(0,1,0.05)
  graph.vars <- sapply(tseq, function(x) calcGraphFeatures(weightThresholding(adj, x)$adj))
  return(graph.vars)
}

weightThresholding <- function(adj, w=0.5){
  # Apply weight-based thresholding to adjacency matrix
  adj[adj < w] <- 0

  return(list(adj=adj, thresh=w))
} 

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
