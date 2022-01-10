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

genIndivNetwork <- function (n, p, q, alpha, obeta0, delta) {
  # n=sparams$n; p=sparams$p; q=sparams$q; zeta=sparams$zeta; alpha=sparams$alpha; delta=sparams$delta;obeta0=sparams$obeta0; beta0=sparams$beta0;xbeta=sparams$xbeta; gbeta = sparams$gbeta; tseq=thresh.seq
  
  
  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, GN, GE for each subject i separately: GN: graph nodes, GE: graph edges
  i=1; X=NULL; GN=NULL; GE=NULL
  
  # ---- (1) Network Generation
  repeat {
    
    # Covariates X_i
    xi=MASS::mvrnorm(1,mu=rep(0,q),diag(delta,q))
    
    # Omega_i: Precision matrix, Kappa_i:mean
    alpha.icpt = alpha[[1]]
    alpha.wmat = alpha[[2]]
    ox= alpha.icpt + c(xi%*%alpha.wmat)
    ox=-ox
    kx=c(xi%*%zeta)
    Omegai=VecToSymMatrix(obeta0, -ox)
    
    # No covariance matrix is generated that is singular
    if(!is.positive.definite(Omegai)){Omegai<-make.positive.definite(Omegai, tol=1e-3)}
    
    # Covariance matrix: Inverse of Omegai
    Sigmai=solve(Omegai)
    
    ### mean matrix ???
    #mui=Sigmai%*%kx
    
    ### generate biomarker nodes M
    #gn=MASS::mvrnorm(1, mui, Sigmai)
    gn=MASS::mvrnorm(1, kx, Sigmai)
    
    ## Partial correlation - Network edge weights
    count=0; sr=1
    ge=numeric((p-1)*p/2); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        # rho^2=pho
        pho=-ox[sr]/sqrt(obeta0[s]*obeta0[r])
        if (abs(pho)>0.9999){count=count+1}
        pho=ifelse(abs(pho)>0.9999, sign(pho)*0.9999, pho)
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
    cc.uw <- clustcoeff(adj, weighted = F)$CC  
    
    graph <- graph_from_adjacency_matrix(adj, diag = F, weighted = T, mode="undirected")
    cpl <- mean_distance(graph, directed = F, unconnected = TRUE)
    cl <- clusters(graph)
    mod <- modularity(graph, membership = cl$membership)
    ass <- assortativity.degree(graph)
    dia <- diameter(graph)
    rad <- radius(graph)
    ev <- eigen_centrality(graph, scale=T, weights=NULL)
    eigen.score <- mean(ev$vector)

  out <- c("cc.w"=cc.w, "cc.uw"=cc.uw, "cpl"=cpl, "mod"=mod, "ass"=ass, "dia"=dia, "rad"=rad, "es"=eigen.score)
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

weightThresholding <- function(adj, w=0.5){
  # Apply weight-based thresholding to adjacency matrix
  adj[adj < w] <- 0

  return(list(adj=adj, thresh=w))
} 

evalIndivNetwork <- function(eweights, msize, tseq, toMatrix=T){
  # eweights=data.graph$GE[1,]; msize=p; tseq=thresh.seq; toMatrix=T
  
  # Perform weight-based thresholding and network feature computation (wrapper function)
  if(toMatrix){adj <- VecToSymMatrix(diag.entry = 0, side.entries = unlist(eweights), mat.size = msize)
  }else{adj <- eweights}
  
  graph.vars <- data.frame(sapply(tseq, function(x){
    out.w <- calcGraphFeatures(weightThresholding(adj, x)$adj)
    out.d <- calcGraphFeatures(densityThresholding(adj, x)$adj)
    out <- c(out.w,out.d)
  } )) %>%
    `colnames<-`(c(tseq)) %>%
    mutate(ThreshMethod=rep(c("weight-based", "density-based"), each=8),
           "Gvar"=str_remove(rownames(.), ".1")) %>%
    melt(id.vars=c("ThreshMethod", "Gvar")) %>%
    `colnames<-`(c("ThreshMethod", "variable", "Thresh", "value")) %>%
    mutate_at(3:4, to_numeric)
  return(graph.vars)
}



# 
#   # ------ (2) Network features
#   #### Sparsification ###
#   MI = abs(MI)
# 
#   # CC for threshold sequence
#   list.gvars <- lapply(1:nrow(MI), function(x) evalSSN(eweights=MI[x,], msize=p, tseq=tseq))
#   data.gvars <- do.call(rbind, lapply(1:length(list.gvars), function(x){
#       colnames(list.gvars[[x]]) <- tseq
#       df.x <- list.gvars[[x]] %>%
#         data.frame() %>%
#         mutate(ThreshMethod=rep(c("weight-based", "density-based"), each=8),
#                "Gvar"=rownames(.)) %>%
#         melt(id.vars=c("ThreshMethod", "Gvar")) %>%
#         `colnames<-`(c("ThreshMethod", "variable", "Thresh", "value")) %>%
#         mutate(Thresh= str_remove(Thresh, "X")) %>%
#         mutate(Subj=x)
#       return(df.x)
#     } )) 
#   
#   df <- list("FMI"=data.gvars, "MI"=MI, "X"=X)
#   
#   return(df)   
# }

