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

to_numeric <- function(x){
  return(as.numeric(as.character(x)))
}

round_0 <- function(x, digits){
  return(sprintf(paste0('%.',digits,'f'),x))
}

# ------------------- MODELLING -----------------------------

get_error_fitted = function(yhat, y) {
  mean.hat <- apply(yhat,1, function(x) mean(x, na.rm = T))
  bias = mean.hat - y
  relbias = ((mean.hat - y)/y)*100
  var = apply(yhat,1, function(x) var(x,na.rm=T))
  
  rmse = apply(((yhat - y) ^ 2), 1, function(x) sqrt(mean(x, na.rm = T)))
  
  out <- cbind(bias, relbias, rmse, var, "mean"=mean.hat)
  return(out)
}


get_error_coef = function(xhat, x) {
  mean.hat <- mean(xhat, na.rm = T)
  bias = mean.hat - x
  relbias = ((mean.hat - x)/x)*100
  var = var(xhat,na.rm=T)
  
  rmse = sqrt(mean((xhat - x) ^ 2), na.rm=T)
  
  out <- c(bias, relbias)
  return(out)
}


# ------------------- NETWORK -----------------------------

genIndivNetwork <- function (n, p, q, omega, mu, delta) {
  # n=sparams$n; p=sparams$p; q=sparams$q; mu=sparams$mu; omega=sparams$pmega; delta=sparams$delta
  
  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, GN, GE for each subject i separately: GN: graph nodes, GE: graph edges
  i=1; X=NULL; GN=NULL; GE=NULL
  
  # ---- (1) Network Generation
  repeat {
    
    # Covariates X_i
    xi=MASS::mvrnorm(1,mu=rep(0,q),diag(delta,q))
    
    # Omega_i: Precision matrix, Kappa_i:mean
    omega.icpt = omega[[1]]
    omega.wmat = omega[[2]]
    ox= omega.icpt + c(xi%*%omega.wmat)
    obeta0=rep(max(abs(ox)),p)   
    
    Mui=c(xi%*%mu)
    Omegai=VecToSymMatrix(obeta0, -ox)
    
    # No covariance matrix is generated that is singular
    if(!is.positive.definite(Omegai)){Omegai<-make.positive.definite(Omegai, tol=1e-3)}
    
    # Covariance matrix: Inverse of Omegai
    Sigmai=solve(Omegai)
    
    ### mean matrix ???
    #mui=Sigmai%*%Mui
    
    ### generate biomarker nodes M
    gn=MASS::mvrnorm(1, Mui, Sigmai)
    
    ## Partial correlation - Network edge weights
    count=0; sr=1
    ge=numeric((p-1)*p/2); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        # rho^2=pho
        pho=-ox[sr]/sqrt(obeta0[s]*obeta0[r])
        if (abs(pho)>0.9999){count=count+1}
        ge[sr]=pho
        sr=sr+1
      }
    }
    
    X=rbind(X,xi)      # covariates
    GN=rbind(GN,gn)    # graph nodes
    GE=rbind(GE,ge)   # graph edges
    
    if (i==n) break
    i=i+1
  }
  
  return(list(X=X, GN=GN, GE=GE))
}


calcGraphFeatures <- function(adj, weighted=NULL){

    cc.w <- clustcoeff(abs(adj), weighted = T)$CC # Clustering coefficient
    graph <- graph_from_adjacency_matrix(adj, diag = F, weighted = T, mode="undirected")
    cpl <- mean_distance(graph, directed = F, unconnected = TRUE)
    cl <- clusters(graph)
    mod <- modularity(graph, membership = cl$membership)
    ass <- assortativity.degree(graph)
    dia <- diameter(graph)
    # ev <- eigen_centrality(graph, scale=T, weights=NULL)
    # eigen.score <- mean(ev$vector)
    ncon <- sum((adj[row(adj)!=col(adj)]) != 0)/2

  out <- c("cc.w"=cc.w, "cpl"=cpl, "mod"=mod, "ass"=ass, "dia"=dia, "ncon"=ncon)
  return(out)
}


densityThresholding <- function(adj, d=0.5, method="trim"){
  # Apply density-based thresholding to adjacency matrix
  V= ncol(adj)
  E.vals=sort(adj[upper.tri(adj)],decreasing = T)
  E.poss <- (V*(V-1))/2
  E.d <- round(E.poss*d,0)
  min.ew <- min(E.vals[1:E.d])
  
  if(method=="trim"){
    repl <- adj[adj> min.ew & row(adj)!=col(adj)]
  }else if(method=="bin"){
    repl <- 1
  }else if(method=="resh"){
    repl <- adj[adj> min.ew & row(adj)!=col(adj)]-min((adj[adj> min.ew & row(adj)!=col(adj)]))
  }else{
    stop("No valid weight replacement method (bin, trim, resh) selected.")
  }
  
  adj[adj <= min.ew] <- 0
  adj[adj > min.ew & row(adj)!=col(adj)] <- repl
  
  return(list(adj=as.matrix(adj), thresh=d))
} 

weightThresholding <- function(adj, w=0.5, method="trim"){
  # Apply weight-based thresholding to adjacency matrix
  # adj=mnet
  
  if(method=="trim"){
    repl <- adj[adj> w & row(adj)!=col(adj)]
  }else if(method=="bin"){
    repl <- 1
  }else if(method=="resh"){
    repl <- adj[adj> w & row(adj)!=col(adj)]-w
  }else{
    stop("No valid weight replacement method (bin, trim, resh) selected.")
  }
  adj[adj <= w] <- 0
  adj[adj > w & row(adj)!=col(adj)] <- repl

  return(list(adj=adj, thresh=w))
} 

wrapperThresholding <- function(eweights, msize, tseq, toMatrix=T){
  # eweights=data.graph$GE[1,]; msize=p; tseq=thresh.seq; toMatrix=T
  # eweights=mnet; msize=p; tseq=thresh.seq; toMatrix=F
  
  # Perform weight-based thresholding and network feature computation (wrapper function)
  if(toMatrix){adj <- VecToSymMatrix(diag.entry = 0, side.entries = unlist(eweights), mat.size = msize)
  }else{adj <- eweights}
  
  thresh.meths = c("trim", "bin", "resh")
  list.vars <- lapply(tseq, function(x){
    out.w <- sapply(thresh.meths, function(y) calcGraphFeatures(weightThresholding(adj, w=x, method = y)$adj)) %>% 
      melt() %>%
      mutate("SparsMethod"="weight-based",
             "Thresh"=x) 
    out.d <- sapply(thresh.meths, function(y) calcGraphFeatures(densityThresholding(adj, d=x, method = y)$adj)) %>% 
      melt() %>%
      mutate("SparsMethod"="density-based",
             "Thresh"=x) 
    out <- data.frame(rbind(out.w,out.d)) %>%
      `colnames<-`(c("Variable", "ThreshMethod", "Value", "SparsMethod", "Thresh")) %>%
      mutate_at(c(3,5), to_numeric)
    return(out)
  } )
  df.vars <- data.frame(do.call(rbind, list.vars))
  return(df.vars)
}


