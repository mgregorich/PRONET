# ============================================================================ #
# Author: MG
# Date: 18.10.2021
# Info: Helper functions for the simulation study
# ============================================================================ #



# ============================ GENERAL =========================================

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

rowwise_scaling01 <- function(x, eps=0.01, ...){
  tmp <- x
  if(length(tmp)>3){
    minx <- ifelse(min(tmp)-eps<0, min(tmp)-eps, 0)
    maxx <- ifelse(max(tmp)>1, max(tmp), 1)
    
    tmp_scaled <- scale(tmp, center = minx, scale = maxx - minx)
    return(tmp_scaled)
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
  return(sqrt(mean((obs-pred)^2, na.rm=T)))
}

calc_rsq <- function(obs, pred){
  return(ifelse(sd(pred)!=0, suppressWarnings(cor(obs,pred)^2), 0))
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
  Brier <- mean((pred-obs)^2, na.rm=T)
  return(Brier)
}

calc_deviance <- function(obs, pred){
  Deviance <- -2*sum( obs*log(pred) + (1-obs)*log(1-pred) )
  return(Deviance)
}

get_error_fitted <- function(yhat, y) {
  mean.hat <- apply(yhat,1, function(x) mean(x, na.rm = T))
  bias <- mean.hat-y
  relbias <- ((mean.hat-y)/y)*100
  var <- apply(yhat,1, function(x) var(x,na.rm=T))
  
  rmse <- apply(((yhat-y)^2), 1, function(x) sqrt(mean(x, na.rm = T)))
  
  out <- cbind(bias, relbias, rmse, var, "mean"=mean.hat)
  return(out)
}


get_error_coef <- function(xhat, x) {
  mean.hat <- mean(xhat, na.rm = T)
  bias <- mean.hat-x
  relbias <- ((mean.hat-x)/x)*100
  var <- var(xhat, na.rm=T)
  
  rmse <- sqrt(mean((xhat-x)^2), na.rm=T)
  
  out <- c(bias, relbias)
  return(out)
}

evalLM <- function(df, k=5, adjust=FALSE){
  # Perform univariable linear regression with CV
  # df=data.AVG$data[[2]]; k=5
  # df=data.oracle$data[[1]]; k=5
  # df=data.null$data[[1]]; k=5
  
  df$fitted <- NA
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("Thresh", "RMSE", "R2", "CS")
  model.form <- as.formula(ifelse(adjust, "Y~Value+X", "Y~Value"))
  
  for(i in 1:k){
    df.train <- df[df$fold !=i, ]
    df.test <- df[df$fold ==i, ]
    
    fit.lm <- lm(model.form, data = df.train, na.action = "na.exclude")
    df.test$fitted <- suppressWarnings(predict(fit.lm, newdata=df.test))
    
    inner[i,] <- c("Thresh"=NA, 
                   "RMSE"=calc_rmse(df.test$Y, df.test$fitted),
                   "R2"=calc_rsq(df.test$Y, df.test$fitted),
                   "CS"=calc_cs(df.test$Y, df.test$fitted))
  }
  return(data.frame("adjust"=adjust,"Thresh"=NA, "RMSE"=mean(inner$RMSE), "R2"=mean(inner$R2), "CS"=mean(inner$CS)))
}


evalLM_dCV <- function(df, k=5, adjust=FALSE){
  # Perform univariable linear regression with double CV for threshold selection 
  # df=data.OPT$data[[1]]; k=5

  df$fitted <- NA
  tseq <- sort(unique(df$Thresh))
  outer_opt <- matrix(NA, nrow=k, ncol=4)
  model.form <- as.formula(ifelse(adjust, "Y~Value+X", "Y~Value"))
  
  for(i in 1:k){
    df.train.out <- df[df$fold !=i, ]
    df.test.out <- df[df$fold ==i, ]
    for(j in unique(df.train.out$fold)){
      inner_opt <- matrix(NA, ncol=2, nrow=length(tseq))
        for(l in 1:length(tseq)){
          df.train.in <- df.train.out[df.train.out$fold !=j & df.train.out$Thresh == tseq[l], ]
          df.test.in <- df.train.out[df.train.out$fold ==j & df.train.out$Thresh == tseq[l], ]
          
          fit.lm <- lm(model.form, data = df.train.in, na.action = "na.exclude")
          df.test.in$fitted <- suppressWarnings(predict(fit.lm, newdata=df.test.in))
          inner_opt[l,] <- c("Thresh"=tseq[l], "RMSE"=calc_rmse(df.test.in$Y, df.test.in$fitted))
      }
    }
      outer_opt[i,1] <- inner_opt[which.min(inner_opt[,2]),1]
      df.train.opt <- df[df$fold !=i & df$Thresh == outer_opt[i,1],]
      df.test.opt <- df[df$fold ==i & df$Thresh == outer_opt[i,1],]
      
      fit.lm <- lm(model.form, data = df.train.opt, na.action = "na.exclude")
      df.test.opt$fitted <- suppressWarnings(predict(fit.lm, newdata=df.test.opt))
      outer_opt[i,2:4] <- c("RMSE"=calc_rmse(df.test.opt$Y, df.test.opt$fitted),
                               "R2"=calc_rsq(df.test.opt$Y, df.test.opt$fitted), 
                               "CS"=calc_cs(df.test.opt$Y, df.test.opt$fitted))
  }
  Thresh = as.numeric(names(sort(table(outer_opt[,1]), decreasing = T))[1])
  RMSE = mean(outer_opt[,2], na.rm=T)
  R2 = mean(outer_opt[,3], na.rm=T)
  CS = mean(outer_opt[,4], na.rm=T)
  return(data.frame("adjust"=adjust,"Thresh"=Thresh,"RMSE"=RMSE, "R2"=R2, "CS"=CS))
}

evalPFR <- function(df, k=5, adjust=FALSE, bs.type="ps", nodes=25, fx=F){
  # Perform scalar-on-function regression with CV
  # df=data.FLEX$data[[22]]; k=5; bs.type="ps"; nodes=NULL; adjust=F; fx=F
  
  df$fitted <- NA
  tmp <- as.matrix.data.frame(df[,str_starts(colnames(df), pattern = "T_")])
  df$M <- tmp[,colSums(is.na(tmp)) < nrow(tmp)/4]
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("Thresh", "RMSE", "R2", "CS")
  tseq <- as.numeric(str_remove(colnames(tmp), "T_"))
  model.form <- as.formula(ifelse(adjust, "Y ~ X + lf(M, k=nodes, bs=bs.type, fx=F)", "Y ~ lf(M, k=nodes, bs=bs.type, fx=fx)"))
  outer_nodes <- matrix(rep(NA, k*2), nrow=5)
  for (i in 1:k){
    df.train.out <- df[df$fold !=i, ]
    df.test.out <- df[df$fold ==i, ]

    fit.fda <- pfr_new(model.form, data=df.train.out, 
                       family="gaussian", method = "REML")
    df.test.out$fitted = c(predict(fit.fda, newdata=df.test.out, type="response"))   
    
    inner[i,] <- c("Thresh"=NA, 
                   "RMSE"=calc_rmse(df.test.out$Y, df.test.out$fitted),
                   "R2"=calc_rsq(df.test.out$Y, df.test.out$fitted),
                   "CS"=calc_cs(df.test.out$Y, df.test.out$fitted))
  } 
  fit.main <- pfr_new(model.form, data=df, 
                      family="gaussian", method="REML")
  
  sd.Xt <- apply(tmp,2, function(x) sd(x, na.rm=T))
  fit.loess <- loess(sd.Xt ~tseq, na.action = "na.exclude")
  sd.Xt.pred <- predict(fit.loess, newdata = coef(fit.main)$X.argvals)
  
  out <- tibble("adjust"=adjust,"Thresh"=NA,"RMSE"=mean(inner$RMSE), "R2"=mean(inner$R2), "CS"=mean(inner$CS), 
                "Coef"=coef(fit.main), "sd.Xt.pred"=sd.Xt.pred) %>%
    nest(Coef=!c(adjust, Thresh, RMSE, R2, CS))
  return(out)
}

# ================================ NETWORK =====================================


genDefaultNetwork <- function(p, q, network.model, beta.params, alpha0.params, alpha12.params, Z1.params, Z2.params){
  
  if(network.model== "scale-free"){
    # Barabasi-Albert model with linear preferential attachment; density > 75% !
    n_edges = 20
    edens = 0
    while(edens < .95){
      BA.graph <- sample_pa(n=p, power=1, m=n_edges, directed = F)
      edens <- edge_density(BA.graph)
      n_edges = n_edges + 1
    }
  }else if(network.model=="small-world"){
    nei_par = p/3
    edens = 0
    while(edens < .95){
      BA.graph <- sample_smallworld(dim=1, size=p, nei=nei_par, p=0.5)
      edens <- edge_density(BA.graph)
      nei_par = nei_par + 1
      }
  }else{
    BA.graph <- sample_gnp(n=p, p=0.95)
    edens <- edge_density(BA.graph)
  }

  BA.strc <- as.matrix(as_adjacency_matrix(BA.graph))

  # -- Edge weights ~ beta(a,b)
  po = (p-1)*p/2                                                                 
  eta.params <- calc_eta_mean_and_var(alpha0.params=alpha0.params, 
                                      Z1.params=Z1.params, Z2.params=Z2.params,
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
  mu = matrix(0,q, p)
  sweight = seq(-2.5,2.5,0.5)
  mu[,sample(1:p, round(p*0.6))] <- sample(sweight, round(p*0.6)*q, replace = T)
  
  return(list("alpha"=alpha, "mu"=mu, "eta.params"=eta.params, "BA.graph"=BA.graph))
}

calc_eta_mean_and_var <- function(alpha0.params, alpha12.params, Z1.params, Z2.params){
  #' alpha0~N(mean1,std1), alpha1~U(a,b), alpha2~U(a.b)
  #' Z1~N(mean2,sd2), Z2~B(1,p)
  #'  Compute mean and std from a linear combination of uniformly and normally distributed variables
  
  alpha12.unif.mean = (alpha12.params[2]-alpha12.params[1])/2
  alpha12.unif.var = ((1/12)*(alpha12.params[2]-alpha12.params[1])^2)
  
  V = alpha0.params[2]^2 + ((alpha12.unif.var+alpha12.unif.mean^2) * (Z1.params[2]^2+Z1.params[1]^2))-
                          (alpha12.unif.mean^2*Z2.params[1]^2) + (alpha12.unif.var+alpha12.unif.mean^2) -
  alpha12.unif.mean^2 * Z2.params[1]
  #S_noisy = V + eps.g^2
  S = sqrt(V)
  M = alpha0.params[1] + alpha12.unif.mean * Z1.params[1] + alpha12.unif.mean * Z2.params[1]
  
  out <- list("mean"=as.numeric(M), "std"=as.numeric(S))
  return(out)
}

transform_to_beta <- function(eta, beta_pars, eta_pars){
  # eta=etai; beta.pars = distr.params$beta; eta.pars = eta.params
  # Convert normally distributed random variable to beta distribution
  p = pnorm(eta, mean=eta_pars[1], sd=eta_pars[2])
  q = qbeta(p, beta_pars[1], beta_pars[2])
  
  # isZero = rbinom(n=length(p), size=1, prob=.75)
  # q <- ifelse(isZero==1, 0, q)


  return(q)
}


genIndivNetwork <- function (n, p, q, eps.g, alpha, Z1.params, Z2.params, mu, beta.params, eta.params) {
  #' Generate n pxp ISN based on BA graph altered by q latent processes

  ## number of possible undirected edges
  po = (p-1)*p/2
  
  ###  Generate X, GN, GE for each subject i separately: GN: graph nodes, GE: graph edges
  eps <- rnorm(n, mean=0, sd=eps.g)
  
  # --- Latent processes z_i
  Z1 <- rnorm(n, mean=Z1.params[1], sd=Z1.params[2])
  Z2 <- rbinom(n, size=1, prob=Z2.params)
  Z = cbind(Z1, Z2)
  
  # --- Edge weights
  eta = alpha[[1]] + outer(alpha[[2]][1,], Z1) + outer(alpha[[2]][2,], Z2)
  teta = eta
  teta[teta!=0] = transform_to_beta(eta=eta[eta!=0], beta_pars = beta.params, eta_pars = unlist(eta.params))

  if(eps.g != 0){
    eta.err = t(t(eta) + eps)
    teta.err = eta.err
    teta.err[teta.err!=0] = transform_to_beta(eta=eta.err[eta.err!=0], beta_pars = beta.params, eta_pars = unlist(eta.params))
  }else{
    teta.err <- teta
  }

  #obeta0 = rep(1,p) 
  #levelplot(teta[1:150,1:150]- teta.err[1:150,1:150])

  
  # ---- (1) Network Generation
  # for(i in 1:n){
    
    #--- Generate biomarker node variables
    #     Mui=c(zi%*%mu)
    # Omegai <- cpp_vec_to_mat(tetai, p)
    
    # Positive definite matrix?
    # if(!is.positive.definite(Omegai)){Omegai2<-make.positive.definite(Omegai, tol=1e-3)}
    
    # Covariance matrix: Inverse of Omegai
    # Sigmai=solve(Omegai)
    
    # Mean matrix
    # mui=Sigmai%*%Mui
    
    # Nodes M
    # gn=MASS::mvrnorm(1, Mui, Sigmai)
    
    #--- Partial correlation - Network edge weights
    # sr =1; ge=numeric((p-1)*p/2); ge.err=numeric((p-1)*p/2);
    # for (s in 1:(p-1)) {
    #   for (r in (s+1):p) {
    #     ge[sr]=tetai[sr]/(obeta0[s]*obeta0[r])
    #     ge.err[sr]=tetai.err[sr]/(obeta0[s]*obeta0[r])
    #     sr=sr+1
    #   }
    # }
    
    # GE[i,] = ge    # graph edges
    # GE.err[i,] = ge.err    # graph edges
    # GN[i,] = gn    # graph nodes
  #}
  
  return(list(Z=Z, GN=0, GE=t(teta), GE.err=t(teta.err)))
}


calcGraphFeatures_new <- function(vec, msize){
  # vec=wt_mat[1,]; msize =p

  cc.w = 1
  cc.uw <- cpp_cc_func(vec, p=msize)
  return(c("cc"=cc.w ,"cc.uw"=cc.uw))
}

calcGraphFeatures <- function(vec, msize){
  
  adj <- VecToSymMatrix(diag.entry = 0, side.entries = vec, mat.size = msize)
  adj[adj>0] <- 1
  
 # cc.w=1
  cc.uw <- mean(WGCNA::clusterCoef(adjMat = adj))
  return(cc.uw)
}

cc_func <- function(A, weighted=F){

    if(!weighted){A[A>0] <- 1
    }else{
      maxh1 = max(as.dist(A))
      minh1 = min(as.dist(A))
      if (maxh1 > 1 | minh1 < 0) 
        stop(paste("The adjacency matrix contains entries that are larger than 1 or smaller than 0: max =", 
                   maxh1, ", min =", minh1))
    }
    diag(A) = 0
    n = ncol(A)
    nolinksNeighbors <- c(rep(-666, n))
    total.edge <- c(rep(-666, n))
    CC <- rep(-666, n)
    
    nolinksNeighbors <- apply(A, 1, function(x) crossprod(x,A) %*% x)
    plainsum <- rowSums(A)
    squaresum <- rowSums(A^2)
    total.edge = plainsum^2 - squaresum
    CCi = ifelse(total.edge == 0, 0, nolinksNeighbors/total.edge)
    CC = mean(CCi)
    return(CC)
}


wrapperThresholding <- function(df, msize, step.size){
  # df=data.network; msize=p; step.size = step.size
  tseq <- seq(0, 1, step.size)

  cc <- cpp_wrapper_thresholding(as.matrix(df), p=msize, step_size=step.size)
  cc <- do.call(rbind, cc)
  res <- data.frame("ID"=rep(1:nrow(df), times=length(tseq)*2),
                    "SparsMethod"=rep(rep(c("weight-based", "density-based"), each=nrow(df)), length(tseq)), 
                    "ThreshMethod"="trim",
                    "Thresh"=rep(tseq, each=nrow(df)*2), 
                    "Variable"="cc.uw",
                    "Value"=cc)
  
  return(res)
}

Thresholding <- function(mat, w=0.5, method="trim", density=F){
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
    if(w!=0){mat[mat >= w] <- 1}
    return(mat)
    
  }else if(method=="resh"){
    mat[mat < w] <- 0
    mat <- t(apply(mat, 1, function(x) rowwise_scaling01(x)))
    return(mat)
    
  }else{
    stop("Select only bin, trim or resh!")
  }
} 

# ===================== REFUND pfr() modification =============================

pfr_new <- function (formula = NULL, fitter = NA, method = "REML", ...) 
{
  require(mgcv)
  
  call <- match.call()
  dots <- list(...)
  if (length(dots)) {
    validDots <- if (!is.na(fitter) && fitter == "gamm4") {
      c(names(formals(gamm4)), names(formals(lmer)))
    }
    else {
      c(names(formals(gam)), names(formals(gam.fit)))
    }
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if (length(notUsed)) 
      warning("Arguments <", paste(notUsed, collapse = ", "), 
              "> supplied but not used.")
  }
  tf <- terms.formula(formula, specials = c("s", "te", "t2", 
                                            "lf", "af", "lf.vd", "re", "peer", "fpc"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]], 
                  simplify = FALSE)
  frmlenv <- environment(formula)
  specials <- attr(tf, "specials")
  where.af <- specials$af - 1
  where.lf <- specials$lf - 1
  where.pr <- specials$peer - 1
  where.fp <- specials$fpc - 1
  where.s <- specials$s - 1
  where.te <- specials$te - 1
  where.t2 <- specials$t2 - 1
  where.re <- specials$re - 1
  where.lf.vd <- specials$lf.vd - 1
  where.all <- c(where.af, where.lf, where.s, where.te, where.t2, 
                 where.re, where.lf.vd, where.pr, where.fp)
  if (length(trmstrings)) {
    where.par <- which(!(1:length(trmstrings) %in% where.all))
  }
  else where.par <- numeric(0)
  responsename <- attr(tf, "variables")[2][[1]]
  newfrml <- paste(responsename, "~", sep = "")
  newfrmlenv <- new.env()
  evalenv <- if ("data" %in% names(call)) 
    eval.parent(call$data)
  else NULL
  nobs <- length(eval(responsename, envir = evalenv, enclos = frmlenv))
  if (missing(fitter) || is.na(fitter)) {
    fitter <- ifelse(nobs > 1e+05, "bam", "gam")
  }
  fitter <- as.symbol(fitter)
  if (as.character(fitter) == "bam" && !("chunk.size" %in% 
                                         names(call))) {
    call$chunk.size <- max(nobs/5, 10000)
  }
  if (as.character(fitter) == "gamm4") 
    stopifnot(length(where.te) < 1)
  assign(x = deparse(responsename), value = as.vector(t(eval(responsename, 
                                                             envir = evalenv, enclos = frmlenv))), envir = newfrmlenv)
  newtrmstrings <- attr(tf, "term.labels")
  if (!attr(tf, "intercept")) {
    newfrml <- paste(newfrml, "0", sep = "")
  }
  where.refund <- c(where.af, where.lf, where.lf.vd, where.pr, 
                    where.fp, where.re)
  if (length(where.refund)) {
    fterms <- lapply(terms[where.refund], function(x) {
      eval(x, envir = evalenv, enclos = frmlenv)
    })
    newtrmstrings[where.refund] <- sapply(fterms, function(x) {
      safeDeparse(x$call)
    })
    lapply(fterms, function(x) {
      lapply(names(x$data), function(nm) {
        assign(x = nm, value = x$data[[nm]], envir = newfrmlenv)
        invisible(NULL)
      })
      if ("xt" %in% names(x$call)) {
        xtvars <- all.vars(x$call$xt)
        if (length(xtvars)) {
          sapply(xtvars, function(xtvar) {
            xtvarval <- eval(as.name(xtvar), envir = evalenv, 
                             enclos = frmlenv)
            assign(x = xtvar, value = xtvarval, envir = parent.env(newfrmlenv))
            invisible(NULL)
          })
        }
      }
      invisible(NULL)
    })
    fterms <- lapply(fterms, function(x) x[names(x) != "data"])
  }
  else fterms <- NULL
  where.mgcv <- c(where.par, where.s, where.te, where.t2)
  if (length(where.mgcv)) {
    if ("data" %in% names(call)) 
      frmlenv <- list2env(eval.parent(call$data), frmlenv)
    lapply(terms[where.mgcv], function(x) {
      nms <- if (!is.null(names(x))) {
        all.vars(x[names(x) == ""])
      }
      else all.vars(x)
      sapply(nms, function(nm) {
        stopifnot(length(get(nm, envir = frmlenv)) == 
                    nobs)
        assign(x = nm, value = get(nm, envir = frmlenv), 
               envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
  }
  newfrml <- formula(paste(newfrml, paste(newtrmstrings, collapse = "+")))
  environment(newfrml) <- newfrmlenv
  pfrdata <- list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv), function(x) {
    if (is.numeric(x) | is.logical(x)) {
      mean(x)
    }
    else if (is.factor(x)) {
      names(which.max(table(x)))
    }
    else NA
  }, simplify = FALSE)
  newcall <- expand.call(pfr, call)
  newcall$fitter <- newcall$bs.int <- newcall$bs.yindex <- NULL
  newcall$formula <- newfrml
  newcall$data <- quote(pfrdata)
  newcall$method <- method
  newcall[[1]] <- fitter
  res <- eval(newcall)
  res.smooth <- if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth
  }
  else res$smooth
  names(res.smooth) <- sapply(res.smooth, function(x) x$label)
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$smooth <- res.smooth
  }
  else {
    res$smooth <- res.smooth
  }
  termtype <- rep("par", length(terms))
  for (i in 1:length(specials)) termtype[specials[[i]] - 1] <- names(specials)[i]
  ret <- list(formula = formula, responsename = responsename, 
              nobs = nobs, termnames = names(terms), termtype = termtype, 
              datameans = datameans, ft = fterms)
  if (as.character(fitter) %in% c("gamm4", "gamm")) {
    res$gam$pfr <- ret
    class(res$gam) <- c("pfr", class(res$gam))
  }
  else {
    res$pfr <- ret
    class(res) <- c("pfr", class(res))
  }
  return(res)
}


safeDeparse <- function(expr){
  # turn an expression into a _single_ string, regardless of the expression's length
  ret <- paste(deparse(expr), collapse="")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}


list2df <- function(l){
  # make a list into a dataframe -- matrices are left as matrices!
  nrows <- sapply(l, function(x) nrow(as.matrix(x)))
  stopifnot(length(unique(nrows)) == 1)
  ret <- data.frame(rep(NA, nrows[1]))
  for(i in 1:length(l)) ret[[i]] <- l[[i]]
  names(ret) <- names(l)
  return(ret)
}

expand.call <- function(definition=NULL, call=sys.call(sys.parent(1)),
                        expand.dots = TRUE)
{
  call <- match.call(definition, call, expand.dots)
  #given args:
  ans <- as.list(call)
  this_function <- safeDeparse(ans[[1]])
  # rm namespace operators so that function finding works (mlr-org/mlr/#1559)
  if (grepl(":", this_function)) {
    this_function <- gsub("(^[^:]+:+)(.+$)", "\\2", this_function)
    if (!exists(this_function)) {
      warning("expand.call couldn't find ", this_function,
              " and returned unchanged call:\n <", call, ">")
      return(call)
    }
  }
  frmls <- formals(this_function)
  #remove formal args with no presets:
  frmls <- frmls[!sapply(frmls, is.symbol)]
  add <- which(!(names(frmls) %in% names(ans)))
  return(as.call(c(ans, frmls[add])))
}



