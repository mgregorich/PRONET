# ============================================================================ #
# Author: MG
# Date: 18.10.2021
# Info: Helper functions for the simulation study
# ============================================================================ #



# ============================ GENERAL =========================================

cbind_results <- function(x, sim.path){
  list.tmp <- readRDS(here::here(sim.path, x))
  tmp <- data.frame(cbind(list.tmp$scenario, list.tmp$tbl_results))
  return(tmp)
  }

restricted_rnorm <- function(n, mean = 0, sd = 1, min = 0, max = 1) {
  # Generalized restricted normal
  bounds <- pnorm(c(min, max), mean, sd)
  u <- runif(n, bounds[1], bounds[2])
  q <- qnorm(u, mean, sd)
  return(q)
}

scaling01 <- function(x, ...){
  # Scale vector between [0,1]
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

# ============================ MODELLING =======================================

c_index <- function(pred, obs){
  c.model <- concreg(data=data.frame(predicted=pred, observed=obs), observed~predicted, npar=TRUE)
  return(1-cindex(c.model))
}
calc_rmse <- function(obs, pred){
  return(sqrt(mean((obs-pred)^2)))
}
calc_rsq <- function(obs, pred){
  return(suppressWarnings(cor(obs,pred)^2))
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

evalLM <- function(data.lm, k=5){
  # Perform univariable linear regression with CV
  # data.lm=data.AVG$data[[2]]; k=5
  # data.lm=data.oracle$data; k=5
  
  df=data.frame(Y=data.lm$Y, X=data.lm$Value, fold=data.lm$fold, fitted=NA)
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("Thresh", "RMSE", "R2", "CS")
  
  for(i in 1:k){
    df.train <- df[df$fold !=i, ]
    df.test <- df[df$fold ==i, ]
    
    fit.lm <- lm(Y~X, data = df.train, na.action = "na.exclude")
    df.test$fitted <- suppressWarnings(predict(fit.lm, newdata=df.test))
    
    inner[i,] <- c("Thresh"=NA, 
                   "RMSE"=calc_rmse(df.test$Y, df.test$fitted),
                   "R2"=calc_rsq(df.test$Y, df.test$fitted),
                   "CS"=calc_cs(df.test$Y, df.test$fitted))
  }
  return(tibble("Thresh"=NA,"RMSE"=mean(inner$RMSE), "R2"=mean(inner$R2), "CS"=mean(inner$CS)))
}


evalLM_dCV <- function(data.lm, k=5){
  # Perform univariable linear regression with double CV for threshold selection 
  # data.lm=data.bRMSE$data[[1]]; k=5

  df=data.frame(Y=data.lm$Y, X=data.lm$Value, Thresh=data.lm$Thresh, fold=data.lm$fold, fitted=NA)
  obest.thresh <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(obest.thresh) <- c("bThresh", "RMSE", "R2", "CS")

    for(i in 1:k){
      df.train <- df[df$fold !=i, ]
      df.test <- df[df$fold ==i, ]
    
      ibest.thresh <- matrix(NA, nrow=5, ncol = 2)
      for(j in unique(df.train$fold)){
        df.train2 <- df.train[df.train$fold !=j, ]
        df.test2 <- df.train[df.train$fold ==j, ]
        
        int.res <- data.frame(t(sapply(unique(df$Thresh), function(x){
          fit.lm <- lm(Y~X, data=df.train2[df.train2$Thresh==x,], na.action = "na.exclude")
          df.test2[df.test2$Thresh==x,]$fitted <- suppressWarnings(predict(fit.lm, newdata=df.test2[df.test2$Thresh==x,]))
          out <- data.frame("Thresh"=x, "RMSE"=calc_rmse(df.test2[df.test2$Thresh==x,]$Y, df.test2[df.test2$Thresh==x,]$fitted))
          return(out)
        })))
        ibest.thresh[j,] <- unlist(int.res[which.min(int.res$RMSE),])
      }
      obt <- ibest.thresh[which.min(ibest.thresh[,2]),1]
      fit.lm <- lm(Y~X, data = df.train[df.train$Thresh==obt,], na.action = "na.exclude")
      df.test[df.test$Thresh==obt,]$fitted <- suppressWarnings(predict(fit.lm, newdata=df.test[df.test$Thresh==obt,]))
      obest.thresh[i,] <- c("best.threshold"= obt, 
                            "RMSE"=calc_rmse(df.test[df.test$Thresh==obt,]$Y, df.test[df.test$Thresh==obt,]$fitted),
                            "R2"=calc_rsq(df.test[df.test$Thresh==obt,]$Y, df.test[df.test$Thresh==obt,]$fitted),
                            "CS"=calc_cs(df.test[df.test$Thresh==obt,]$Y, df.test[df.test$Thresh==obt,]$fitted))
    }
  
  Thresh = as.numeric(names(sort(table(obest.thresh$bThresh)))[1])
  RMSE = mean(obest.thresh$RMSE, na.rm=T)
  R2 = mean(obest.thresh$R2, na.rm=T)
  CS = mean(obest.thresh$CS, na.rm=T)
  return(tibble("Thresh"=Thresh,"RMSE"=RMSE, "R2"=R2, "CS"=CS))
}

evalPFR <- function(data.fda, k=5, tseq){
  # Perform scalar-on-function regression with CV
  # data.fda=data.FDA$data[[1]]; k=5
  
  flength=length(tseq)
  df=data.frame("fold"=data.fda$fold, "Y"=data.fda$Y, "fitted"=NA)
  df$X <- as.matrix.data.frame(data.fda[,(ncol(data.fda)-flength+1):ncol(data.fda)])
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("Thresh", "RMSE", "R2", "CS")
  
  for (i in 1:k){
    df.train <- df[df$fold !=i, ]
    df.test <- df[df$fold ==i, ]
    
    fit.fda <- refund::pfr(Y ~ lf(X, k = 10, bs="ps"), data=df.train, family="gaussian")
    df.test$fitted = c(predict(fit.fda, newdata=df.test, type="response"))   
    
    inner[i,] <- c("Thresh"=NA, 
                   "RMSE"=calc_rmse(df.test$Y, df.test$fitted),
                   "R2"=calc_rsq(df.test$Y, df.test$fitted),
                   "CS"=calc_cs(df.test$Y, df.test$fitted))
  } 
  
  fit.main <- refund::pfr(Y ~ lf(X, k = 10, bs="ps"), data=df, family="gaussian")

  out <- tibble("Thresh"=NA,"RMSE"=mean(inner$RMSE), "R2"=mean(inner$R2), "CS"=mean(inner$CS), "Coef"=coef(fit.main)) %>%
    nest(Coef=!c(Thresh, RMSE, R2, CS))
  return(out)
}

# ================================ NETWORK =====================================


genDefaultNetwork <- function(p, q, BA.graph, beta.params, alpha0.params, alpha12.params, X.params){
  
  # -- Barabasi-Albert model for Bernoulli graph
  BA.strc <- as.matrix(as_adjacency_matrix(BA.graph))

  # -- Edge weights ~ beta(a,b)
  po = (p-1)*p/2                                                                 
  eta.params <- calc_eta_mean_and_var(alpha0.params=alpha0.params, 
                                      X.params=X.params,
                                      alpha12.params=alpha12.params)
  alpha0 <- BA.strc[lower.tri(BA.strc)]
  alpha0[alpha0==1] <- rnorm(sum(alpha0), alpha0.params[1], alpha0.params[2])
  omega.imat=matrix(alpha0,1, po, byrow = T)
  
  alpha12 <- rep(BA.strc[lower.tri(BA.strc)], q)
  alpha12[alpha12==1] <- runif(sum(alpha12), alpha12.params[1], alpha12.params[2])                             
  alpha12.wmat=matrix(alpha12, q, po, byrow = T)
  alpha=list("alpha0"=alpha0, "alpha12"=alpha12.wmat)
  
  
  # -- Only important in case variables are generated on which the network can be estimated
  # mu: qxp matrix of weighting for mean
  mu=matrix(0,q, p)
  sweight=seq(-2.5,2.5,0.5)
  mu[,sample(1:p, round(p*0.6))] <- sample(sweight, round(p*0.6)*q, replace = T)
  
  return(list("alpha"=alpha, "mu"=mu, "eta.params"=eta.params))
}

calc_eta_mean_and_var <- function(alpha0.params, alpha12.params, X.params){
  #' alpha0~N(mean1,std1), alpha1~U(a,b), alpha2~U(a.b)
  #' X1~N(mean2,sd2), X2~N(mean2,sd2)
  #'  Compute mean and std from a linear combination of uniformly and normally distributed variables
  
  alpha12.unif.mean = (alpha12.params[2]-alpha12.params[1])/2
  alpha12.unif.var = ((1/12)*(alpha12.params[2]-alpha12.params[1])^2)
  
  V=alpha0.params[2]^2 + 2*(((alpha12.unif.var +alpha12.unif.mean^2)*(X.params[2]^2+X.params[1]^2))-
                                    (alpha12.unif.mean^2*X.params[1]^2))
  S=sqrt(V)
  M=alpha0.params[1] + alpha12.unif.mean * X.params[1] + alpha12.unif.mean * X.params[1]
  
  out <- list("mean"=as.numeric(M), "sd"=as.numeric(S))
  return(out)
}

transform_to_beta <- function(eta, beta.pars, eta.pars){
  # eta=etai; beta.pars = distr.params$beta; eta.pars = eta.params
  # Convert normally distributed random variable to beta distribution
  p = pnorm(eta, mean=eta.pars$mean, sd=eta.pars$sd)
  q = qbeta(p, beta.pars[1], beta.pars[2])
  return(q)
}


genIndivNetwork <- function (n, p, q, alpha, X.params, mu, beta.params, eta.params) {
  #' Generate n pxp ISN based on BA graph altered by q latent processes

  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, GN, GE for each subject i separately: GN: graph nodes, GE: graph edges
  i=1; X=NULL; GN=NULL; GE=NULL
  
  # ---- (1) Network Generation
  repeat {
    
    # Covariates X_i
    xi=MASS::mvrnorm(1,mu=rep(X.params[1],q),diag(X.params[2]^2,q))
    
    # Omega_i: Precision matrix, Kappa_i:mean
    alpha0 = alpha[[1]]; alpha12 = alpha[[2]]
    etai = alpha0 + c(xi%*%alpha12)
    ox=etai
    ox[ox!=0] = transform_to_beta(eta=etai[etai!=0], beta.pars = beta.params, eta.pars = eta.params)

    Mui=c(xi%*%mu)
    obeta0 = rep(1,p) 
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
    d <- edge_density(graph)
    # cl <- components(graph)
    # mod <- modularity(graph, membership = cl$membership)
    # ass <- assortativity.degree(graph)
    # dia <- diameter(graph)
    # ev <- eigen_centrality(graph, scale=T, weights=NULL)
    # eigen.score <- mean(ev$vector)
    ncon <- sum((adj[row(adj)!=col(adj)]) != 0)/2

  out <- c("cc"=cc.w, "cpl"=cpl, "ncon"=ncon, "dens"=d)
  return(out)
}


densityThresholding <- function(adj, d=0.5, method="trim"){
  # Apply density-based thresholding to adjacency matrix
 
  adj.new <- adj
  
  E.vals=sort(adj[upper.tri(adj)],decreasing = T)
  E.vals <- E.vals[!E.vals<=0]
  E.d <- round(length(E.vals)*d,0)
  min.ew <- min(E.vals[1:E.d])

  if(method=="trim"){
    repl <- adj[adj>= min.ew & row(adj)!=col(adj)]
  }else if(method=="bin"){
    repl <- 1
  }else if(method=="resh"){
    repl <- adj[adj>= min.ew & row(adj)!=col(adj)]-min(0,adj[adj>= min.ew & row(adj)!=col(adj)])
  }else{
    stop("No valid weight replacement method (bin, trim, resh) selected.")
  }
  
  adj.new[adj < min.ew] <- 0
  adj.new[adj >= min.ew & row(adj)!=col(adj)] <- repl
  
  return(list(adj=as.matrix(adj.new), thresh=d))
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
  # eweights=data.network[2,]; msize=p; tseq=da.thresh; toMatrix=T
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





