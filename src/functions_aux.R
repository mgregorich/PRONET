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
  y <- (x-min(x, ...))/(max(x, ...)-min(x, ...))
  
  return(y)}

rowwise_scaling01 <- function(x, w, ...){
  tmp <- x[x >= w]
  if(length(tmp)>3){
    minx <- min(tmp)
    maxx <- max(tmp)
    
    tmp_scaled <- scale(tmp, center = minx, scale = maxx - minx)
    x[x >= w] <- tmp_scaled
    return(x)
  }else{return(x)}
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

to_numeric <- function(x){
  return(as.numeric(as.character(x)))
}

round_0 <- function(x, digits){
  return(sprintf(paste0('%.',digits,'f'),x))
}

source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}

# ============================ MODELLING =======================================

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

calc_cstat <- function(obs, pred, outcome_typ="gaussian"){
  
  if(outcome_typ=="binomial"){
    obs <- as.numeric(as.character(obs))
    cstat <- somers2(pred, obs)["C"]
    return(cstat)
  }else{
    if(sd(pred)!=0){
      if(cor(pred, obs, use = 'pairwise.complete.obs') > 0.98){
        return(cor(pred, obs, use = 'pairwise.complete.obs'))
      }else{
        c.model <- concreg(data=data.frame(predicted=pred, observed=obs), observed~predicted, npar=TRUE, maxit=200)
        return(1-cindex(c.model))
      }
    }else{return(NA)}    
  }
}

calc_brier <- function(obs, pred){
  obs <- as.numeric(as.character(obs))
  Brier <- mean((pred-obs)^2)
  return(Brier)
}

calc_deviance <- function(obs, pred){
  Deviance<- -2*sum( obs*log(pred) + (1-obs)*log(1-pred) )
  return(Deviance)
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
  # data.lm=data.oracle$data[[1]]; k=5
  # data.lm=data.null$data[[1]]; k=5
  
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
  return(tibble("Thresh"=NA, "RMSE"=mean(inner$RMSE), "R2"=mean(inner$R2), "CS"=mean(inner$CS)))
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
          df.test2[df.test2$Thresh==x,]$fitted <- predict(fit.lm, newdata=df.test2[df.test2$Thresh==x,])
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

evalPFR <- function(data.fda, k=5, bs.type="ps", nodes=20){
  # Perform scalar-on-function regression with CV
  # data.fda=data.FDA$data[[1]]; k=5; bs.type="ps"; nodes=20
  
  df=data.frame("fold"=data.fda$fold, "Y"=as.numeric(as.character(data.fda$Y)), "fitted"=NA)
  tmp <- as.matrix.data.frame(data.fda[,str_starts(colnames(data.fda), pattern = "T_")])
  df$X <- tmp[,colSums(is.na(tmp)) < nrow(tmp)/4]
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("Thresh", "RMSE", "R2", "CS")
  tseq <- as.numeric(str_remove(colnames(tmp), "T_"))
  
  for (i in 1:k){
    df.train <- df[df$fold !=i, ]
    df.test <- df[df$fold ==i, ]
    
    fit.fda <- refund::pfr(Y ~ lf(X, k = nodes, bs=bs.type), data=df.train, family="gaussian")
    df.test$fitted = c(predict(fit.fda, newdata=df.test, type="response"))   
    
    inner[i,] <- c("Thresh"=NA, 
                   "RMSE"=calc_rmse(df.test$Y, df.test$fitted),
                   "R2"=calc_rsq(df.test$Y, df.test$fitted),
                   "CS"=calc_cs(df.test$Y, df.test$fitted))
  } 
  
  fit.main <- refund::pfr(Y ~ lf(X, k = nodes, bs=bs.type), data=df, family="gaussian")

  sd.Xt <- apply(tmp,2, function(x) sd(x, na.rm=T))
  fit.loess <- loess(sd.Xt ~tseq, na.action = "na.exclude")
  sd.Xt.pred <- predict(fit.loess, newdata = coef(fit.main)$X.argvals)
  
  out <- tibble("Thresh"=NA,"RMSE"=mean(inner$RMSE), "R2"=mean(inner$R2), "CS"=mean(inner$CS), 
                "Coef"=coef(fit.main), "sd.Xt"=sd.Xt.pred) %>%
    nest(Coef=!c(Thresh, RMSE, R2, CS))
  return(out)
}


# ================================ NETWORK =====================================


genDefaultNetwork <- function(p, q, BA.graph, beta.params, alpha0.params, alpha12.params, X1.params, X2.params){
  
  # -- Barabasi-Albert model for Bernoulli graph
  BA.strc <- as.matrix(as_adjacency_matrix(BA.graph))

  # -- Edge weights ~ beta(a,b)
  po = (p-1)*p/2                                                                 
  eta.params <- calc_eta_mean_and_var(alpha0.params=alpha0.params, 
                                      X1.params=X1.params, X2.params=X2.params,
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

calc_eta_mean_and_var <- function(alpha0.params, alpha12.params, X1.params, X2.params){
  #' alpha0~N(mean1,std1), alpha1~U(a,b), alpha2~U(a.b)
  #' X1~N(mean2,sd2), X2~B(1,p)
  #'  Compute mean and std from a linear combination of uniformly and normally distributed variables
  
  alpha12.unif.mean = (alpha12.params[2]-alpha12.params[1])/2
  alpha12.unif.var = ((1/12)*(alpha12.params[2]-alpha12.params[1])^2)
  
  V=alpha0.params[2]^2 + ((alpha12.unif.var +alpha12.unif.mean^2)*(X1.params[2]^2+X1.params[1]^2))-
                                    (alpha12.unif.mean^2*X2.params[1]^2) + (alpha12.unif.var+alpha12.unif.mean^2) -
    alpha12.unif.mean^2 * X2.params[1]
  S_noisy = V + eps.g^2
  S=sqrt(V)
  M=alpha0.params[1] + alpha12.unif.mean * X1.params[1] + alpha12.unif.mean * X2.params[1]
  
  out <- list("mean"=as.numeric(M), "std"=as.numeric(S))
  return(out)
}

transform_to_beta <- function(eta, beta.pars, eta.pars){
  # eta=etai; beta.pars = distr.params$beta; eta.pars = eta.params
  # Convert normally distributed random variable to beta distribution
  p = pnorm(eta, mean=eta.pars$mean, sd=eta.pars$std)
  q = qbeta(p, beta.pars[1], beta.pars[2])
  return(q)
}


genIndivNetwork <- function (n, p, q, alpha, X1.params, X2.params, mu, beta.params, eta.params) {
  #' Generate n pxp ISN based on BA graph altered by q latent processes

  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, GN, GE for each subject i separately: GN: graph nodes, GE: graph edges
  i=1; X=NULL; GN=NULL; GE=NULL
  
  # ---- (1) Network Generation
  repeat {
    
    # Covariates X_i
    xi_1=rnorm(q/2, mean=X1.params[1], sd=X1.params[2])
    xi_2=rbinom(q/2, size=1, prob=X2.params)
    xi <- c(xi_1, xi_2)
    
    # Omega_i: Precision matrix, Kappa_i:mean
    alpha0 = alpha[[1]]; alpha12 = alpha[[2]]
    etai = alpha0 + c(xi%*%alpha12)
    ox=etai
    ox[ox!=0] = transform_to_beta(eta=etai[etai!=0], beta.pars = beta.params, eta.pars = eta.params)
   # hist(ox)
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
    # sr =1;ge=numeric((p-1)*p/2);
    # for (s in 1:(p-1)) {
    #   for (r in (s+1):p) {
    #     pho=ox[sr]/(obeta0[s]*obeta0[r])
    #     ge[sr]=log(1-pho)/2
    #     sr=sr+1
    #   }
    # }
    ge =ox
   # hist(ge)
    
    X=rbind(X,xi)      # covariates
    GN=rbind(GN,gn)    # graph nodes
    GE=rbind(GE,ge)   # graph edges
    
    if (i==n) break
    i=i+1
  }
  
  return(list(X=X, GN=GN, GE=GE))
}


calcGraphFeatures <- function(adj, weighted=NULL, w=0){
  
  #swi <- qgraph::smallworldness(adj, B=10)[1]
  cc.w <- mean(WGCNA::clusterCoef(adj))
  cc.uw <- clustcoeff(adj, weighted = F)[[1]]
  
  #graph <- graph_from_adjacency_matrix(adj, diag = F, weighted = T, mode="undirected")
  # cpl <- mean_distance(graph, directed = F, unconnected = TRUE)
  # ass <- assortativity.degree(graph)
  # m<-wsyn::cluseigen(adj)
  # mod<-wsyn::modularity(adj, m[[length(m)]], decomp=F)
  
  # AUC until x
  # d_fun <- ecdf(adj[upper.tri(adj)])
  # auc <- 1-d_fun(w)
  
  # Mean edge weight
  # mean_ew <- mean(adj[upper.tri(adj)])
  
  # edens <- edge_density(graph)
  
  # out <- c("cc"=cc.w , "cpl"=cpl,  "ass"=ass, "mod"=mod, "auc"=auc, "edens"=edens)
  out <- c("cc"=cc.w ,"cc.uw"=cc.uw)
  
  return(out)
}


densityThresholding <- function(adj, d=0.5, method="trim"){
  # Apply density-based thresholding to adjacency matrix
  # method = "resh"; d=0.1; 
  
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
    tmp <- adj[adj>= min.ew & row(adj)!=col(adj)] 
    repl <- ifelse(length(tmp)>3, scaling01(tmp), tmp)
  }else{
    stop("No valid weight replacement method (bin, trim, resh) selected.")
  }
  
  adj.new[adj < min.ew] <- 0
  adj.new[adj >= min.ew & row(adj)!=col(adj)] <- repl
  
  return(list(adj=as.matrix(adj.new), thresh=d))
} 

weightThresholding <- function(adj, w=0.5, method="trim"){
  # Apply weight-based thresholding to adjacency matrix
  # adj=adj; method="trim"; w=0.33
  adj.new <- adj
  
  if(method=="trim"){
    repl <- adj[adj>= w & row(adj)!=col(adj)]
  }else if(method=="bin"){
    repl <- 1
  }else if(method=="resh"){
    tmp <- adj[adj>= w & row(adj)!=col(adj)]
    if(length(tmp)>3){repl <- scaling01(tmp)
    }else{ repl <- tmp}
  }else{
    stop("No valid weight replacement method (bin, trim, resh) selected.")
  }
  adj.new[adj < w] <- 0
  adj.new[adj >= w & row(adj)!=col(adj)] <- repl

  return(list(adj=adj.new, thresh=w))
} 

calc_gvars <- function(vec, msize){
  # vec=wt_mat[1,]; msize=p
  
  adj <- VecToSymMatrix(diag.entry = 0, side.entries = unlist(vec), mat.size = msize)
  cc.w <- mean(WGCNA::clusterCoef(adj))
  cc.uw <- clustcoeff(adj, weighted = F)[[1]]
  
  return(c("cc"=cc.w ,"cc.uw"=cc.uw))
}

wrapper_func <- function(df, msize){
  # df=data.network; msize=p

  tmeth = c("trim", "resh", "bin")
  tseq = seq(0, 1, 0.02)
  
  res <- list()
  i = 1
  for(t in tseq){
    for(m in tmeth){
      wt_mat <- thresh_func(mat=df, w=t, method = m, density = F)
      dt_mat <- thresh_func(mat=df, w=t, method = m, density = T)
      
      res_wt <- t(apply(wt_mat, 1, function(x) c("SparsMethod"="weight-based","ThreshMethod"=m, "Thresh"=t, 
                                                 calc_gvars(x, msize = msize))))
      res_dt <- t(apply(dt_mat, 1, function(x) c("SparsMethod"="density-based","ThreshMethod"=m, "Thresh"=t, 
                                                 calc_gvars(x, msize = msize))))
      
      res_tmp <- data.frame("Subj"=1:nrow(df), rbind(res_wt, res_dt))
      res[[i]] <- melt(res_tmp, id.vars = c("Subj", "SparsMethod", "ThreshMethod", "Thresh"), 
                       variable.name ="Variable", value.name = "Value")
      i = i + 1
    }
  }
  return(res)
}

thresh_func <- function(mat, w=0.5, method="trim", density=F){
  # Apply weight-based thresholding to adjacency matrix
  # mat=df; method="bin"; w=0

  if(density){
    E.vals <- rowSort(as.matrix(mat), descending=T)
    E.d <- round(ncol(E.vals)*w,0)
    w <- ifelse(E.d!=0, unlist(E.vals[,E.d]), 1)
  }
  
  if(method=="trim"){
    mat[mat < w] <- 0
    return(mat)
    
  }else if(method=="bin"){
    mat[mat < w] <- 0
    mat[mat >= w] <- 1
    return(mat)
    
  }else if(method=="resh"){
    mat[mat < w] <- 0
    mat <- t(apply(mat, 1, function(x) rowwise_scaling01(x,w=w)))
    return(mat)
    
  }else{
    stop("Select only bin, trim or resh!")
  }
} 


wrapperThresholding <- function(eweights, msize, toMatrix=T){
  # eweights=data.network[2,]; msize=p; toMatrix=T
  # eweights=mnet; msize=p; toMatrix=F
  
  # Perform weight-based thresholding and network feature computation (wrapper function)
  if(toMatrix){adj <- VecToSymMatrix(diag.entry = 0, side.entries = unlist(eweights), mat.size = msize)
  }else{adj <- eweights}
  
  thresh.meths = c("trim", "bin", "resh")
  tseq = seq(0, 1, 0.02)
  list.vars.w <- lapply(tseq, function(x){
    out.w <- lapply(thresh.meths, function(y) calcGraphFeatures(weightThresholding(adj, w=x, method = y)$adj, w=x)) 
    out.w <- do.call(cbind, out.w, quote=TRUE) %>%
      t() %>%
      data.frame() %>%
      mutate("SparsMethod"="weight-based",
             "ThreshMethod"=thresh.meths,
             "Thresh"=x) %>%
      pivot_longer(cols=-c(SparsMethod, ThreshMethod, Thresh), values_to = "Value", names_to = "Variable")
    return(out.w)
  } )
  list.vars.d <- lapply(tseq, function(x){
    out.d <- lapply(thresh.meths, function(y) calcGraphFeatures(densityThresholding(adj, d=x, method = y)$adj, w=x)) 
    out.d <- do.call(cbind, out.d, quote=TRUE) %>% 
      t() %>%
      data.frame() %>%
      mutate("SparsMethod"="density-based",
             "ThreshMethod"=thresh.meths,
             "Thresh"=x) %>%
      pivot_longer(cols=-c(SparsMethod, ThreshMethod, Thresh), values_to = "Value", names_to = "Variable")
    return(out.d)
  } )      
  list.vars <- c(list.vars.w, list.vars.d)
  df.vars <- as.data.frame(do.call(rbind, list.vars, quote=TRUE))
  return(df.vars)
}




