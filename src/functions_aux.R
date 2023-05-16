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

transform_to_beta <- function(eta, beta_pars){
  # eta=etai; beta.pars = distr.params$beta; eta.pars = eta.params
  # Convert normally distributed random variable to beta distribution
  pemp <- ecdf(eta)
  p = pemp(eta)
  q = qbeta(p, beta_pars[1], beta_pars[2])
  
  return(q)
}

scaling01 <- function(x, uplim=1, ...){
  # Scale vector between [0,1]
  y <- (x-min(x, ...))/((max(x, ...)-min(x, ...)*uplim))
  
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
  # dat=data.oracle$data[[2]]; k=5; family="gaussian"; adjust=T
  
  if(!any(family==c("gaussian", "binomial"))){stop("family must be gaussian or binomial")}
  
  dat$fitted <- NA
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("data_ana_t", "metric_1", "metric_2", "metric_3")
  model.form <- as.formula(ifelse(adjust, "Y~value+X", "Y~value"))
  
  for(i in 1:k){
    dat.train <- dat[dat$fold !=i, ]
    dat.test <- dat[dat$fold ==i, ]
    
    fit.tmp <- glm(model.form, data = dat.train, na.action = "na.exclude", family=family)
    dat.test$fitted <- suppressWarnings(predict(fit.tmp, newdata=dat.test, type="response"))
    dat[dat$fold ==i, ]$fitted <- dat.test$fitted
    
    inner[i,] <- c("data_ana_t"=NA, 
                   "metric_1"=ifelse(family=="gaussian", 
                                     calc_rmspe(dat.test$Y, dat.test$fitted), 
                                     calc_brier(dat.test$Y, dat.test$fitted)),
                   "metric_2"=calc_rsq(dat.test$Y, dat.test$fitted),
                   "metric_3"=ifelse(family=="gaussian", 
                                     calc_cs(dat.test$Y, dat.test$fitted), 
                                     calc_cstat(dat.test$Y, dat.test$fitted, outcome_typ="binomial")))
    
  }
  fit.main <- glm(model.form, data = dat, na.action = "na.exclude", family=family)
  
  out <- tibble("adjust"=adjust,
                "data_ana_t"=NA, 
                "metric_1"=mean(inner$metric_1, na.rm = T), 
                "metric_2"=mean(inner$metric_2, na.rm = T), 
                "metric_3"=mean(inner$metric_3, na.rm = T),
                "coef"=nest(tibble("coef"=coef(fit.main)), data=everything()),
                "pred"=nest(tibble("Y"=dat[!is.na(dat$fitted),]$Y, "fitted"= dat[!is.na(dat$fitted),]$fitted), data=everything())) %>%
    unnest(coef) %>%
    dplyr::rename("coef"=data) %>%
    unnest(pred) %>%
    dplyr::rename("pred"=data) 
  
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
  colnames(obest.thresh) <- c("bThresh", "metric_1", "metric_2", "metric_3")
  model.form <- as.formula(ifelse(adjust, "Y~value+X", "Y~value"))
  
  for(i in 1:k){
    dat.train <- dat[dat$fold !=i, ]
    dat.test <- dat[dat$fold ==i, ]
    
    ibest.thresh <- matrix(NA, nrow=k, ncol = 2)
    for(j in unique(dat.train$fold)){
      dat.train2 <- dat.train[dat.train$fold !=j, ]
      dat.test2 <- dat.train[dat.train$fold ==j, ]
      
      int.res <- matrix(NA, ncol = 2, nrow=length(unique(dat$data_ana_t)))
      for(l in 1:length(unique(dat$data_ana_t))){
        x = sort(unique(dat$data_ana_t))[l]
        if(sd(dat.train2[dat.train2$data_ana_t==x,]$value, na.rm = T)!=0){
          fit.tmp.in <- glm(model.form, data=dat.train2[dat.train2$data_ana_t==x,], na.action = "na.exclude", family=family)
          dat.test2[dat.test2$data_ana_t==x,]$fitted <- predict(fit.tmp.in, newdata=dat.test2[dat.test2$data_ana_t==x,], type="response")
          eval_metric <- ifelse(family=="gaussian", 
                                calc_rmspe(dat.test2[dat.test2$data_ana_t==x,]$Y, dat.test2[dat.test2$data_ana_t==x,]$fitted),
                                calc_brier(dat.test2[dat.test2$data_ana_t==x,]$Y, dat.test2[dat.test2$data_ana_t==x,]$fitted))
          int.res[l,] <- c("data_ana_t"=x, "metric_1"=eval_metric)
        }else{
          int.res[l,] <- c("data_ana_t"=x, "metric_1"=NA)
        }
      }
      ibest.thresh[j,] <- int.res[which.min(int.res[,2]),]
    }
    opt_t <- ibest.thresh[which.min(ibest.thresh[,2]),1]
    fit.tmp.out <- glm(model.form, data = dat.train[dat.train$data_ana_t==opt_t,], na.action = "na.exclude", family=family)
    dat.test[dat.test$data_ana_t==opt_t,]$fitted <- suppressWarnings(predict(fit.tmp.out, newdata=dat.test[dat.test$data_ana_t==opt_t,], type="response"))
    dat[dat$fold ==i & dat$data_ana_t==opt_t, ]$fitted <- dat.test[dat.test$data_ana_t==opt_t,]$fitted
    obest.thresh[i,] <- c("best.threshold"= opt_t, 
                          "metric_1"=ifelse(family=="gaussian",
                                            calc_rmspe(dat.test[dat.test$data_ana_t==opt_t,]$Y, dat.test[dat.test$data_ana_t==opt_t,]$fitted),
                                            calc_brier(dat.test[dat.test$data_ana_t==opt_t,]$Y, dat.test[dat.test$data_ana_t==opt_t,]$fitted)),
                          "metric_2"=calc_rsq(dat.test[dat.test$data_ana_t==opt_t,]$Y, dat.test[dat.test$data_ana_t==opt_t,]$fitted),
                          "metric_3"=ifelse(family=="gaussian",
                                            calc_cs(dat.test[dat.test$data_ana_t==opt_t,]$Y, dat.test[dat.test$data_ana_t==opt_t,]$fitted),
                                            calc_cstat(dat.test[dat.test$data_ana_t==opt_t,]$Y, dat.test[dat.test$data_ana_t==opt_t,]$fitted,outcome_typ = "binomial")))
  }
  # Either select most often chosen threshold or if no threshold most often, then with smallest metric
  opt_t_final <- ifelse(any(table(obest.thresh$bThresh)>2), as.numeric(names(sort(table(obest.thresh$bThresh)))[1]), obest.thresh[which.min(obest.thresh[,2]),1])
  fit.main <- glm(model.form, data = dat[dat$data_ana_t==as.character(opt_t_final),], na.action = "na.exclude", family=family)
  
  out <- tibble("adjust"=adjust,
                "data_ana_t" = opt_t_final,
                "metric_1" = mean(obest.thresh$metric_1, na.rm=T),
                "metric_2" = mean(obest.thresh$metric_2, na.rm=T),
                "metric_3" = mean(obest.thresh$metric_3, na.rm=T),
                "coef"=nest(tibble("coef"=coef(fit.main)), data=everything()),
                "pred"=nest(tibble("Y"=dat[!is.na(dat$fitted),]$Y, "fitted"= dat[!is.na(dat$fitted),]$fitted), data=everything())) %>%
    unnest(coef) %>%
    dplyr::rename("coef"=data) %>%
    unnest(pred) %>%
    dplyr::rename("pred"=data) 
  
  if(family=="gaussian"){
    colnames(out)[3:5] <- c("RMSPE", "R2", "CS")
  }else{
    colnames(out)[3:5] <- c("Brier", "R2", "C")
  }
  
  return(out)
}

perform_FLEX <- function(data.fda, k=5, adjust=FALSE, bs.type="ps", bs.dim=25, family="gaussian"){
  # Perform scalar-on-function regression with CV
  # data.fda=data.FLEX$data[[1]]; k=5; bs.type="ps"; bs.dim=25; adjust=F; family="gaussian"
  
  if(!any(family==c("gaussian", "binomial"))){stop("family must be gaussian or binomial")}
  
  dat <- data.frame("fold"=data.fda$fold, "Y"=as.numeric(as.character(data.fda$Y)), "fitted"=NA, "X2"=data.fda$X)
  dat$X1 <- as.matrix.data.frame(data.fda[,str_starts(colnames(data.fda), pattern = "T_")])
  
 # point.constraint <- switch(data.fda$SM[1], "density-based"=paste0("pc=1"), "weight-based"=paste0("pc=0"))
  if(adjust){
    # @fx ... fixed regression spline fx=TRUE; penalized spline fx=FALSE
    # @pc ... point constraint; forces function to f(x)=0 at x=1
    model.form <- as.formula(paste0("Y ~ X2 + lf(X1, k = bs.dim, bs=bs.type)")) 
  }else{model.form <- as.formula(paste0("Y ~ lf(X1, k = bs.dim, bs=bs.type)"))}
  
  inner <- data.frame(matrix(NA, nrow=k, ncol=4))
  colnames(inner) <- c("data_ana_t", "metric_1", "metric_2", "metric_3")
  for (i in 1:k){
    dat.train.out <- dat[dat$fold !=i, ]
    dat.test.out <- dat[dat$fold ==i, ]
    
    fit.fda <- pfr_new(model.form, data=dat.train.out, 
                       family=family, method = "REML")
    dat.test.out$fitted <- c(predict(fit.fda, newdata=dat.test.out, type="response"))   
    dat[dat$fold ==i, ]$fitted <- dat.test.out$fitted
    
    inner[i,] <- c("data_ana_t"=NA, 
                   "metric_1"=ifelse(family=="gaussian", 
                                     calc_rmspe(dat.test.out$Y, dat.test.out$fitted), 
                                     calc_brier(dat.test.out$Y, dat.test.out$fitted)),
                   "metric_2"=calc_rsq(dat.test.out$Y, dat.test.out$fitted),
                   "metric_3"=ifelse(family=="gaussian", 
                                     calc_cs(dat.test.out$Y, dat.test.out$fitted), 
                                     calc_cstat(dat.test.out$Y, dat.test.out$fitted, outcome_typ="binomial")))
    
  } 
  fit.main <- pfr_new(model.form, data=dat,  family=family, method="REML")
  
  sd.Xt <- apply(dat$X1, 2, function(x) sd(x, na.rm=T))
  out <- tibble("adjust"=adjust,
                "data_ana_t"=NA,
                "metric_1"=mean(inner$metric_1, na.rm=T), 
                "metric_2"=mean(inner$metric_2, na.rm=T), 
                "metric_3"=mean(inner$metric_3, na.rm=T), 
                "coef"=nest(tibble("coef"=coef(fit.main), "sd.Xt"=sd.Xt), data=everything()),
                "pred"=nest(tibble("Y"=dat$Y, "fitted"= dat$fitted), data=everything())) %>%
    unnest(coef) %>%
    dplyr::rename("coef"=data) %>%
    unnest(pred) %>%
    dplyr::rename("pred"=data) 
  
  if(family=="gaussian"){
    colnames(out)[3:5] <- c("RMSPE", "R2", "CS")
  }else{
    colnames(out)[3:5] <- c("Brier", "R2", "C")
  }
  
  return(out)
}

# ================================ NETWORK =====================================

noisecor <- function(cormat, epsilon = .01, eidim=2, scaling=1){
  #' Modified from https://github.com/MarkJBrewer/ICsims/blob/master/R/noisecor.R
  #' @references 
  #' For full details, see
  #' \cite{Hardin, J., Garcia, S. R., and Golan, D. (2013). A method for generating realistic correlation matrices. Annals of Applied Statistics, 7(3):1733-1762.}
  
  ndim <- dim(cormat)[1]
  diag(cormat) <- 1 - epsilon
  eivect <- c( )
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2)) )
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- cormat + (bigE * scaling)
  cor.nz
}

calcGraphFeatures <- function(adj, weighted=NULL){
  
  cc.w <- mean(WGCNA::clusterCoef(adj))
  # graph <- graph_from_adjacency_matrix(adj, diag = F, weighted = T, mode="undirected")
  # cpl <- mean_distance(graph, directed = F, unconnected = TRUE)
  
  cc.w[is.nan(cc.w)] <- 0
  #cpl[is.nan(cpl)] <- 0
  cpl<-0
  out <- c("cc"=cc.w , "cpl"=cpl)
  return(out)
}


wrapperThresholding <- function(dat, set_ids, mdim, step_size, graph_feature){
  # dat=data.network; mdim=p; step_size = step_size; set_ids=unique(df$data$ID)
  tseq <- seq(0, 1, step_size)
  dat <- apply(as.matrix(dat),2, as.numeric)
  val_graph_feature <- cpp_wrapper_thresholding(dat, p=mdim, step_size=step_size, feature=as.character(graph_feature))
  val_graph_feature <- do.call(rbind, val_graph_feature)
  res <- data.frame("ID"=rep(set_ids, times=length(tseq)*2),
                    "data_ana_thresh"=rep(rep(c("weight-based", "density-based"), each=nrow(dat)), length(tseq)), 
                    "data_ana_t"=rep(tseq, each=nrow(dat)*2), 
                    "variable"=graph_feature,
                    "value"=val_graph_feature )
  
  return(res)
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


