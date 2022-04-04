# ============================================================================ #
# Author: MG
# Date: 18.10.2021
# Info: Helper functions for the simulation study
# ============================================================================ #


# ------------------- GENERAL -----------------------------

restricted_rnorm <- function(n, mean = 0, sd = 1, min = 0, max = 1) {
  # Generalized restricted normal
  bounds <- pnorm(c(min, max), mean, sd)
  u <- runif(n, bounds[1], bounds[2])
  q <- qnorm(u, mean, sd)
  return(q)
}

scaling01 <- function(x, ...){
  y <- (x - min(x, ...)) / (max(x, ...) - min(x, ...))
  return(y)}


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

c_index <- function(pred, obs){
  c.model <- concreg(data=data.frame(predicted=pred, observed=obs), observed~predicted, npar=TRUE)
  return(1-cindex(c.model))
}
calc_rmse <- function(obs, pred){
  return(sqrt(mean((obs-pred)^2)))
}
calc_rsq <- function(obs, pred){
  return(cor(obs,pred)^2)
}
calc_cs <- function(obs,pred){
  if(all(is.na(pred))){
    cs<-NA
  }else{
    cs<-lm(obs~pred, na.action = "na.exclude")$coefficients[2]
  }
  return(cs)
}

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


evalLM <- function(Y, X, fold, k=5){
  # Perform univariable linear regression and extract coeffs, rmse and adjusted r2
  # X <- unlist(data.gvars %>% filter(SparsMethod == "density-based"& ThreshMethod == "trim" & Variable == "cc" & Thresh == 0) %>% select(Value))
  # Y <- data.gvars %>% filter(SparsMethod == "density-based"& ThreshMethod == "trim" & Variable == "cc" & Thresh == 0) %>% select(Y)
  # fold <- data.gvars %>% filter(SparsMethod == "density-based"& ThreshMethod == "trim" & Variable == "cc" & Thresh == 0) %>% select(fold)
  df=data.frame(Y=Y, X=X, fold=fold, fitted=NA)

    for(i in 1:k){
    df.train <- df[df$fold !=i, ]
    df.test <- df[df$fold ==i, ]

    if(!(all(is.na(df$X)) | var(df$X)==0)){
      fit.lm <- lm(Y~X, data = df.train, na.action = "na.exclude")
      df[df$fold==i,]$fitted <- suppressWarnings(predict(fit.lm, newdata=df.test))
      }
    }
  out <- c(to_numeric(df$fitted))
  return(out)
}

evalPFR <- function(x, fold, k=5, tseq){
  # x = data.gvars %>%
  #   pivot_wider(values_from = Value, names_from = Thresh) %>%
  #   filter(SparsMethod == "density-based"& ThreshMethod == "trim" & Variable == "cc")
  # y= data.gvars[!duplicated(data.gvars$Y),]$Y
  # fold= data.gvars[!duplicated(data.gvars$Y),]$fold
  # k=5
  
  flength=length(tseq)
  df=data.frame("fold"=x$fold, "Y"=x$Y, "fitted"=NA)
  df$X <- as.matrix.data.frame(x[,(ncol(x)-flength+1):ncol(x)])
  
  for (i in 1:k){
    df.train <- df[df$fold !=i, ]
    df.test <- df[df$fold ==i, ]
    
    fit.fda <- refund::pfr(Y ~ lf(X, k = 5, bs="ps"), data=df.train, family="gaussian")
    df[df$fold==i,]$fitted = c(predict(fit.fda, newdata=df.test, type="response"))    
  } 
  
  fit.main <- refund::pfr(Y ~ lf(X, k = 5, bs="ps"), data=df, family="gaussian")

  rmse <- calc_rmse(df$Y,df$fitted)
  r2 <- calc_rsq(df$Y, df$fitted)
  cs <- calc_cs(df$Y, df$fitted)
  
  out <- list()
  out$res <- data.frame("RMSE"=rmse, "R2"=r2, "CS"=cs)
  out$coef <- coef(fit.main)
  return(out)
}

# ------------------- NETWORK -----------------------------

calc_eta_mean_and_var <- function(norm.pars, unif.pars){
  # alpha0~N(mean,std), alpha1~U(a,b), alpha2~U(a.b)
  # X1~N(0,1), X2~N(0,1)
  unif.mean = (unif.pars["max"]-unif.pars["min"])/2
  unif.var = ((1/12)*(unif.pars["max"]-unif.pars["min"])^2)
  S=sqrt(norm.pars["sd"]^2 + 2*((unif.var +(unif.mean^2)*1)-(unif.mean^2*0)))
  M=norm.pars["mean"] + unif.mean * 0 + unif.mean * 0
  return(list("mean"=as.numeric(M), "sd"=as.numeric(S)))
}

transform_to_beta <- function(eta, beta.pars, eta.pars){
  # eta=etai; beta.pars = distr.params$beta; eta.pars = eta.params

  p = pnorm(eta, mean=eta.pars$mean, sd=eta.pars$sd)
  q = qbeta(p, beta.pars["shape1"], beta.pars["shape2"])
  return(q)
}


genIndivNetwork <- function (n, p, q, alpha, distr.params, eta.params, mu, delta) {
  # n=sparams$n; p=sparams$p; q=sparams$q; mu=sparams$mu; alpha=sparams$alpha; delta=sparams$delta
  
  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, GN, GE for each subject i separately: GN: graph nodes, GE: graph edges
  i=1; X=NULL; GN=NULL; GE=NULL
  
  # ---- (1) Network Generation
  repeat {
    
    # Covariates X_i
    xi=MASS::mvrnorm(1,mu=rep(distr.params$X.norm[1],q),diag(distr.params$X.norm[2]^2,q))
    
    # Omega_i: Precision matrix, Kappa_i:mean
    alpha0 = alpha[[1]]; alpha12 = alpha[[2]]
    etai = alpha0 + c(xi%*%alpha12)
    ox=etai
    ox[ox!=0] = transform_to_beta(eta=etai[etai!=0], beta.pars = distr.params$beta, eta.pars = eta.params)
    obeta0 = rep(1,p) 

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
    sr=1
    ge=numeric((p-1)*p/2); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        # rho^2=pho
        pho=-ox[sr]/sqrt(obeta0[s]*obeta0[r])
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
    # cl <- components(graph)
    # mod <- modularity(graph, membership = cl$membership)
    # ass <- assortativity.degree(graph)
    # dia <- diameter(graph)
    # ev <- eigen_centrality(graph, scale=T, weights=NULL)
    # eigen.score <- mean(ev$vector)
    ncon <- sum((adj[row(adj)!=col(adj)]) != 0)/2

  out <- c("cc"=cc.w, "cpl"=cpl, "ncon"=ncon)
  return(out)
}


densityThresholding <- function(adj, d=0.5, method="trim"){
  # Apply density-based thresholding to adjacency matrix

  E.vals=sort(adj[upper.tri(adj)],decreasing = T)
  E.vals <- E.vals[!E.vals==0]
  E.d <- round(length(E.vals)*d,0)
  min.ew <- min(E.vals[1:E.d])

  if(method=="trim"){
    repl <- adj[adj> min.ew & row(adj)!=col(adj)]
  }else if(method=="bin"){
    repl <- 1
  }else if(method=="resh"){
    repl <- adj[adj> min.ew & row(adj)!=col(adj)]-min(0,adj[adj> min.ew & row(adj)!=col(adj)])
  }else{
    stop("No valid weight replacement method (bin, trim, resh) selected.")
  }
  
  adj[adj <= min.ew] <- 0
  adj[adj > min.ew & row(adj)!=col(adj)] <- repl
  
  return(list(adj=as.matrix(adj), thresh=d))
} 

weightThresholding <- function(adj, w=0.5, method="trim"){
  # Apply weight-based thresholding to adjacency matrix
  # adj=mat; method="trim"; w=0.33
  
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
  # eweights=data.network[2,]; msize=p; tseq=thresh.seq; toMatrix=T
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





