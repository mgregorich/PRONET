# ============================================================================ #
# Author: MG
# Date: 18.10.2021
# Info: Helper functions for the simulation study
# ============================================================================ #



# ============================ GENERAL =========================================

to_factor <- function(x){
  as.character(as.factor(x))
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

rowwise_scaling01 <- function(x, eps=0.01, ...){
  tmp <- x
  if(length(tmp)>3){
    minx <- ifelse(min(tmp)-eps<0, min(tmp)-eps, 0)
    maxx <- ifelse(max(tmp)>1, max(tmp), 1)
    
    tmp_scaled <- scale(tmp, center = minx, scale = maxx - minx)
    return(tmp_scaled)
  }else{return(x)}
}

VecToSymMatrix <- function(diag.entries, side.entries, mat.size, byrow=T){
  # Generate symmetric matrix from vectors holding diagonal values and side entries
  # diag.entry=0; side.entries=grow; mat.size=150
  side.entries <- unlist(side.entries)
  diags = rep(diag.entries, mat.size)
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

calc_rmspe <- function(obs, pred){
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
  
  RMSPE <- apply(((yhat-y)^2), 1, function(x) sqrt(mean(x, na.rm = T)))
  
  out <- cbind(bias, relbias, RMSPE, var, "mean"=mean.hat)
  return(out)
}


get_error_coef <- function(xhat, x) {
  mean.hat <- mean(xhat, na.rm = T)
  bias <- mean.hat-x
  relbias <- ((mean.hat-x)/x)*100
  var <- var(xhat, na.rm=T)
  
  RMSPE <- sqrt(mean((xhat-x)^2), na.rm=T)
  
  out <- c(bias, relbias)
  return(out)
}


perform_AVG <- function(dat, k=5, adjust=T, family="gaussian"){
  # Perform univariable linear regression with CV
  # dat=data.AVG$data[[3]]; k=3; family="binomial"; adjust=F
  
  if(!any(family==c("gaussian", "binomial"))){stop("family must be gaussian or binomial")}
  
  dat$fitted <- NA
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("Thresh", "Metric_1", "Metric_2", "Metric_3")
  model.form <- as.formula(ifelse(adjust, "Y~Value+X", "Y~Value"))
  
  for(i in 1:k){
    dat.train <- dat[dat$fold !=i, ]
    dat.test <- dat[dat$fold ==i, ]
    
    fit.tmp <- glm(model.form, data = dat.train, na.action = "na.exclude", family=family)
    dat.test$fitted <- suppressWarnings(predict(fit.tmp, newdata=dat.test, type="response"))
    dat[dat$fold ==i, ]$fitted <- dat.test$fitted
    
    inner[i,] <- c("Thresh"=NA, 
                   "Metric_1"=ifelse(family=="gaussian", 
                                     calc_rmspe(dat.test$Y, dat.test$fitted), 
                                     calc_brier(dat.test$Y, dat.test$fitted)),
                   "Metric_2"=calc_rsq(dat.test$Y, dat.test$fitted),
                   "Metric_3"=ifelse(family=="gaussian", 
                                     calc_cs(dat.test$Y, dat.test$fitted), 
                                     calc_cstat(dat.test$Y, dat.test$fitted, outcome_typ="binomial")))
    
  }
  fit.main <- glm(model.form, data = dat, na.action = "na.exclude", family=family)
  
  out <- tibble("Adjust"=adjust,
                "Thresh"=NA, 
                "Metric_1"=mean(inner$Metric_1, na.rm = T), 
                "Metric_2"=mean(inner$Metric_2, na.rm = T), 
                "Metric_3"=mean(inner$Metric_3, na.rm = T),
                "Coef"=nest(tibble("Coef"=coef(fit.main)), data=everything()),
                "Pred"=nest(tibble("Y"=dat[!is.na(dat$fitted),]$Y, "fitted"= dat[!is.na(dat$fitted),]$fitted), data=everything())) %>%
    unnest(Coef) %>%
    dplyr::rename("Coef"=data) %>%
    unnest(Pred) %>%
    dplyr::rename("Pred"=data) 
  
  if(family=="gaussian"){
    colnames(out)[3:5] <- c("RMSPE", "R2", "CS")
  }else{
    colnames(out)[3:5] <- c("Brier", "R2", "C")
  }
  
  return(out)
}


perform_OPT <- function(dat, k=5, adjust=F, family="gaussian"){
  # Perform univariable linear regression with double CV for threshold selection 
  # dat=data.OPT$data[[2]]; k=k; adjust=F; family="gaussian"
  if(!any(family==c("gaussian", "binomial"))){stop("family must be gaussian or binomial")}
  
  dat$fitted <- NA
  #dat$Y <- ifelse(family=="gaussian", as.numeric(dat$Y), as.factor(as.character(dat$Y)))
  obest.thresh <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(obest.thresh) <- c("bThresh", "Metric_1", "Metric_2", "Metric_3")
  model.form <- as.formula(ifelse(adjust, "Y~Value+X", "Y~Value"))
  
  for(i in 1:k){
    dat.train <- dat[dat$fold !=i, ]
    dat.test <- dat[dat$fold ==i, ]
    
    ibest.thresh <- matrix(NA, nrow=k, ncol = 2)
    for(j in unique(dat.train$fold)){
      dat.train2 <- dat.train[dat.train$fold !=j, ]
      dat.test2 <- dat.train[dat.train$fold ==j, ]
      
      int.res <- matrix(NA, ncol = 2, nrow=length(unique(dat$Thresh)))
      for(l in 1:length(unique(dat$Thresh))){
        x = sort(unique(dat$Thresh))[l]
        if(sd(dat.train2[dat.train2$Thresh==x,]$Value, na.rm = T)!=0){
          fit.tmp.in <- glm(model.form, data=dat.train2[dat.train2$Thresh==x,], na.action = "na.exclude", family=family)
          dat.test2[dat.test2$Thresh==x,]$fitted <- predict(fit.tmp.in, newdata=dat.test2[dat.test2$Thresh==x,], type="response")
          eval_metric <- ifelse(family=="gaussian", 
                                calc_rmspe(dat.test2[dat.test2$Thresh==x,]$Y, dat.test2[dat.test2$Thresh==x,]$fitted),
                                calc_brier(dat.test2[dat.test2$Thresh==x,]$Y, dat.test2[dat.test2$Thresh==x,]$fitted))
          int.res[l,] <- c("Thresh"=x, "Metric_1"=eval_metric)
        }else{
          int.res[l,] <- c("Thresh"=x, "Metric_1"=NA)
        }
      }
      ibest.thresh[j,] <- int.res[which.min(int.res[,2]),]
    }
    opt_t <- ibest.thresh[which.min(ibest.thresh[,2]),1]
    fit.tmp.out <- glm(model.form, data = dat.train[dat.train$Thresh==opt_t,], na.action = "na.exclude", family=family)
    dat.test[dat.test$Thresh==opt_t,]$fitted <- suppressWarnings(predict(fit.tmp.out, newdata=dat.test[dat.test$Thresh==opt_t,], type="response"))
    dat[dat$fold ==i & dat$Thresh==opt_t, ]$fitted <- dat.test[dat.test$Thresh==opt_t,]$fitted
    obest.thresh[i,] <- c("best.threshold"= opt_t, 
                          "Metric_1"=ifelse(family=="gaussian",
                                            calc_rmspe(dat.test[dat.test$Thresh==opt_t,]$Y, dat.test[dat.test$Thresh==opt_t,]$fitted),
                                            calc_brier(dat.test[dat.test$Thresh==opt_t,]$Y, dat.test[dat.test$Thresh==opt_t,]$fitted)),
                          "Metric_2"=calc_rsq(dat.test[dat.test$Thresh==opt_t,]$Y, dat.test[dat.test$Thresh==opt_t,]$fitted),
                          "Metric_3"=ifelse(family=="gaussian",
                                            calc_cs(dat.test[dat.test$Thresh==opt_t,]$Y, dat.test[dat.test$Thresh==opt_t,]$fitted),
                                            calc_cstat(dat.test[dat.test$Thresh==opt_t,]$Y, dat.test[dat.test$Thresh==opt_t,]$fitted,outcome_typ = "binomial")))
  }
  # Either select most often chosen threshold or if no threshold most often, then with smallest metric
  opt_t_final <- ifelse(any(table(obest.thresh$bThresh)>2), as.numeric(names(sort(table(obest.thresh$bThresh)))[1]), obest.thresh[which.min(obest.thresh[,2]),1])
  fit.main <- glm(model.form, data = dat[dat$Thresh==as.character(opt_t_final),], na.action = "na.exclude", family=family)
  
  out <- tibble("Adjust"=adjust,
                "Thresh" = opt_t_final,
                "Metric_1" = mean(obest.thresh$Metric_1, na.rm=T),
                "Metric_2" = mean(obest.thresh$Metric_2, na.rm=T),
                "Metric_3" = mean(obest.thresh$Metric_3, na.rm=T),
                "Coef"=nest(tibble("Coef"=coef(fit.main)), data=everything()),
                "Pred"=nest(tibble("Y"=dat[!is.na(dat$fitted),]$Y, "fitted"= dat[!is.na(dat$fitted),]$fitted), data=everything())) %>%
    unnest(Coef) %>%
    dplyr::rename("Coef"=data) %>%
    unnest(Pred) %>%
    dplyr::rename("Pred"=data) 
  
  if(family=="gaussian"){
    colnames(out)[3:5] <- c("RMSPE", "R2", "CS")
  }else{
    colnames(out)[3:5] <- c("Brier", "R2", "C")
  }
  
  return(out)
}

perform_FLEX <- function(data.fda, k=5, adjust=FALSE, bs.type="ps", bs.dim=15, family="gaussian"){
  # Perform scalar-on-function regression with CV
  # data.fda=data.FLEX$data[[1]]; k=5; bs.type="ps"; bs.dim=25; adjust=F; family="gaussian"
  
  if(!any(family==c("gaussian", "binomial"))){stop("family must be gaussian or binomial")}
  
  dat=data.frame("fold"=data.fda$fold, "Y"=as.numeric(as.character(data.fda$Y)), "fitted"=NA, "X2"=data.fda$X)
  tmp <- as.matrix.data.frame(data.fda[,str_starts(colnames(data.fda), pattern = "T_")])
  # dat$X1<- tmp[,colSums(is.na(tmp)) < nrow(tmp)/6]
  dat$X1 <- tmp
  
  point.constraint <- switch(data.fda$SM[1], "density-based"=paste0("pc=0"), "weight-based"=paste0("pc=1"))
  if(adjust){
    # @fx ... fixed regression spline fx=TRUE; penalized spline fx=FALSE
    # @pc ... point constraint; forces function to f(x)=0 at x=1
    model.form <- as.formula(paste0("Y ~ X2 + lf(X1, k = bs.dim, bs=bs.type, ",point.constraint,")")) 
  }else{model.form <- as.formula(paste0("Y ~ lf(X1, k = bs.dim, bs=bs.type, ",point.constraint,")"))}
  
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("Thresh", "Metric_1", "Metric_2", "Metric_3")
  for (i in 1:k){
    dat.train.out <- dat[dat$fold !=i, ]
    dat.test.out <- dat[dat$fold ==i, ]
    
    fit.fda <- pfr_new(model.form, data=dat.train.out, 
                       family=family, method = "REML")
    dat.test.out$fitted <- c(predict(fit.fda, newdata=dat.test.out, type="response"))   
    dat[dat$fold ==i, ]$fitted <- dat.test.out$fitted
    
    inner[i,] <- c("Thresh"=NA, 
                   "Metric_1"=ifelse(family=="gaussian", 
                                     calc_rmspe(dat.test.out$Y, dat.test.out$fitted), 
                                     calc_brier(dat.test.out$Y, dat.test.out$fitted)),
                   "Metric_2"=calc_rsq(dat.test.out$Y, dat.test.out$fitted),
                   "Metric_3"=ifelse(family=="gaussian", 
                                     calc_cs(dat.test.out$Y, dat.test.out$fitted), 
                                     calc_cstat(dat.test.out$Y, dat.test.out$fitted, outcome_typ="binomial")))
    
  } 
  fit.main <- pfr_new(model.form, data=dat,  family=family, method="REML")
  
  sd.Xt <- apply(tmp,2, function(x) sd(x, na.rm=T))
  out <- tibble("Adjust"=adjust,
                "Thresh"=NA,
                "Metric_1"=mean(inner$Metric_1, na.rm=T), 
                "Metric_2"=mean(inner$Metric_2, na.rm=T), 
                "Metric_3"=mean(inner$Metric_3, na.rm=T), 
                "Coef"=nest(tibble("Coef"=coef(fit.main), "sd.Xt"=sd.Xt), data=everything()),
                "Pred"=nest(tibble("Y"=dat$Y, "fitted"= dat$fitted), data=everything())) %>%
    unnest(Coef) %>%
    dplyr::rename("Coef"=data) %>%
    unnest(Pred) %>%
    dplyr::rename("Pred"=data) 
  
  if(family=="gaussian"){
    colnames(out)[3:5] <- c("RMSPE", "R2", "CS")
  }else{
    colnames(out)[3:5] <- c("Brier", "R2", "C")
  }
  
  return(out)
}

# ================================ NETWORK =====================================



genDefaultNetwork <- function(n, p, q, beta.params, alpha0.params, alpha12.params, Z1.params, Z2.params){
  
  po = ((p-1)*p)/2                                                                 
  eta.params <- calc_eta_mean_and_var(alpha0.params=alpha0.params, 
                                      Z1.params=Z1.params, Z2.params=Z2.params,
                                      alpha12.params=alpha12.params)
  
  alpha0.imat <- matrix(NA, ncol=po, nrow = n)
  iedens = runif(n, min=0.75, max=1)
  for(i in 1:n){
    # Barabasi-Albert model with linear preferential attachment; density > 75% !
    n_edges = 20
    edens = 0
    while(edens < iedens[i]){
      default.graph <- sample_pa(n=p, power=1, m=n_edges, directed = F)
      edens <- edge_density(default.graph)
      n_edges = n_edges + 1
    }
    default.strc <- as.matrix(as_adjacency_matrix(default.graph))
    alpha0 <- default.strc[lower.tri(default.strc)]
    
    # -- Edge weights ~ beta(a,b)
    len_alpha0_1 <- sum(alpha0>0)
    initial_weights <- rnorm(len_alpha0_1, alpha0.params[1], alpha0.params[2])
    alpha0[alpha0>0] <- initial_weights
    alpha0.imat[i,]=matrix(alpha0,1, po, byrow = T)    
  }

  alpha12 <- runif(2*length(alpha0), alpha12.params[1], alpha12.params[2])
  alpha12.wmat <- matrix(alpha12, q, po, byrow = T)
  alpha=list("alpha0"=alpha0.imat, "alpha12"=alpha12.wmat)

  return(list("alpha"=alpha, "eta.params"=eta.params, "indv.edensity"=iedens))
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
  
  # V = ((alpha12.unif.var+alpha12.unif.mean^2) * (Z1.params[2]^2+Z1.params[1]^2))-
  #   (alpha12.unif.mean^2*Z2.params[1]^2) + (alpha12.unif.var+alpha12.unif.mean^2) -
  #   alpha12.unif.mean^2 * Z2.params[1]
  # #S_noisy = V + eps.g^2
  # S = sqrt(V)
  # M = alpha12.unif.mean * Z1.params[1] + alpha12.unif.mean * Z2.params[1]
  
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


genIndivNetwork <- function (n, p, q, eps.g, alpha0.params, alpha12.params, Z1.params, Z2.params, beta.params) {
  #' Generate n pxp ISN based on default graph altered by q latent processes

  # -- Setup default network
  dnw.params <- genDefaultNetwork(n, p, q, beta.params, alpha0.params=alpha0.params, alpha12.params=alpha12.params, Z1.params, Z2.params)
  alpha <- dnw.params$alpha
  
  ## number of possible undirected edges
  po = (p-1)*p/2
  eta.params = unlist(dnw.params$eta.params)
  
  ###  Generate X, GN, GE for each subject i separately: GN: graph nodes, GE: graph edges
  # --- Latent processes z_i
  Z1 <- rnorm(n, mean=Z1.params[1], sd=Z1.params[2])
  Z2 <- rbinom(n, size=1, prob=Z2.params)
  Z <- cbind(Z1, Z2)
  
  eta <- alpha[[1]] + t(outer(alpha[[2]][1,], Z1) + outer(alpha[[2]][2,], Z2))
  eta[alpha[[1]]==0] <- 0
  teta <- eta
  teta[teta>0] <- transform_to_beta(eta=teta[teta>0], beta_pars = beta.params, eta_pars = eta.params)

  if(eps.g > 0){
    # eps.tmp <- rnorm(n, mean=0, sd=eps.g)
    # eps <- matrix(rnorm(po, mean=eps.tmp, sd=0), nrow=n, ncol=po)
    eps <- rnorm(n, mean=0, sd=eps.g)
    
    eta.err = eta + eps
    eta.err[alpha[[1]]==0] <- 0
    teta.err = eta.err
    teta.err[eta.err>0] = transform_to_beta(eta=eta.err[eta.err>0], beta_pars = beta.params, eta_pars = eta.params)
  }else{
    eta.err = NA
    teta.err <- teta}
  
  return(list(Z=Z, GN=0, GE=teta, GE.err=teta.err, eta=eta, eta.err=eta.err, indv.edensity=dnw.params$indv.edensity))
}


calcGraphFeatures <- function(vec, msize){
  
  adj <- VecToSymMatrix(diag.entries = 0, side.entries = vec, mat.size = msize)
  adj[adj>0] <- 1
  
 # cc.w=1
  cc.uw <- mean(WGCNA::clusterCoef(adjMat = adj))
  return(cc.uw)
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

sine_fn <- function(x){return(-cos(2*pi*x/0.75)-3/2*sin(2*pi*x/0.75)-2*cos(2*2*pi*x/0.75)+1/2*sin(2*2*pi*x/0.75)+3)/2}

# ===================== REFUND pfr() modification =============================

pfr_new <- function (formula = NULL, fitter = NA, method = "REML", family=NULL, ...) 
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



