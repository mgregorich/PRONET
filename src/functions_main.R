# ============================================================================ #
# Author: MG
# Date: 19.04.2022
# Info: Main simulation functions 
# ============================================================================ #



# ================================== 00. Simulate scenario  =====================================

simulate_scenario <- function(scn){
  #' Given specific design parameters, performs a number of iterations and saves the result in a R object
  # scn = scenarios[19,]
  sourceCpp(here::here("src","utils.cpp"))
  
  file_dgspars = ifelse(scn$dg.spars %in% "weight-based", "wb", "db")
  filename <- paste0("sim_i", scn$iter,"_",scn$outcome, "_",scn$setting, "_", scn$network.model, "_n",scn$n,"_p",scn$p, "_DGS",file_dgspars,"_DGF",names(scn$dg.thresh), 
                     "_beta",unlist(scn$beta.params)[1], "", unlist(scn$beta.params)[2],
                     "_1b",scn$b1,"_2b",scn$b2,"_epsY",scn$epslevel.y,"_epsG",scn$epslevel.g)

  # Preprocess
  beta.params = unlist(scn$beta.params, use.names = F)
  alpha0.params = unlist(scn$alpha0.params, use.names = F)
  alpha12.params = unlist(scn$alpha12.params, use.names = F)
  Z1.params = unlist(scn$Z1.params, use.names = F)
  Z2.params = unlist(scn$Z2.params, use.names = F)

  # -- Setup default network
  dnw.params <- genDefaultNetwork(scn$p, scn$q, scn$network.model, beta.params, alpha0.params, alpha12.params, Z1.params, Z2.params)
  
  # -- Data generation & analysis
  results.sim <- list()
  results.sim <- lapply(1:scn$iter, function(x){
    data.iter <- generate_data(setting = scn$setting,
                               outcome = scn$outcome,
                               n = scn$n, 
                               p = scn$p, 
                               q = scn$q,
                               alpha = dnw.params$alpha, 
                               mu = dnw.params$mu, 
                               eta.params = dnw.params$eta.params,
                               beta.params = beta.params,
                               Z1.params = Z1.params,
                               Z2.params = Z2.params,
                               b0 = scn$b0,
                               b1 = scn$b1,  
                               b2 = scn$b2,
                               eps.y = scn$eps.y, 
                               eps.g = scn$eps.g, 
                               dg.thresh = scn$dg.thresh,
                               dg.spars = scn$dg.spars,
                               step.size = scn$step.size)
    results.iter <- analyse_data(setting = scn$setting,
                                 df = data.iter, 
                                 n = scn$n, 
                                 p = scn$p,
                                 b1 = scn$b1,
                                 dg.thresh = scn$dg.thresh, 
                                 dg.spars = scn$dg.spars,
                                 network.model = scn$network.model,
                                 k = 5,
                                 step.size=scn$step.size) 
  
    return(results.iter)
  })

  
  
    # -- Summarize & save results
  summarize_data(results.sim, setting=scn$setting, outcome=scn$outcome, n=scn$n, p=scn$p, q=scn$q, network.model=scn$network.model,
                 alpha0.params = alpha0.params, alpha12.params = alpha12.params, 
                 Z1.params=Z1.params, Z2.params=Z2.params,beta.params=beta.params, eta.params=eta.params, 
                 b0=scn$b0, b1=scn$b1, b2=scn$b2, eps.y=scn$eps.y, eps.g=scn$eps.g, epslevel.y=scn$epslevel.y, epslevel.g=scn$epslevel.g, 
                 dg.thresh=scn$dg.thresh, dg.spars=scn$dg.spars, default.graph = dnw.params$default.graph, 
                 filename=filename, excel=scn$excel)
}

# ============================ 01. Data generation =============================
generate_data <- function (setting, outcome, n, p, q, mu, alpha, Z1.params, Z2.params, beta.params, eta.params, 
                           b0, b1, b2, eps.y, eps.g, dg.thresh, dg.spars, step.size) {
  #' Data generation (code adapted and modified; initially from https://github.com/shanghongxie/Covariate-adjusted-network)
  # setting=scn$setting; outcome=scn$outcome; n=scn$n; p=scn$p; q=scn$q;
  # alpha=dnw.params$alpha; mu=dnw.params$mu; eta.params = dnw.params$eta.params;
  # beta.params = unlist(scn$beta.params); Z1.params = unlist(scn$Z1.params); Z2.params = unlist(scn$Z2.params);
  # b0=scn$b0; b1 = scn$b1; b2 = scn$b2;
  # eps.y=scn$eps.y; eps.g=scn$eps.g; dg.thresh=scn$dg.thresh; dg.spars=scn$dg.spars

  # -- Individual-specific networks generation
  # Generate ISNs
  po = (p-1)*p/2    
  dg.method <- ifelse(names(dg.thresh) %in% c("flat", "half-sine", "sine"), "func", names(dg.thresh) )
  dg.thresh <- unlist(dg.thresh)
  data.graph <- genIndivNetwork(n=n, p=p, q=q, eps.g=eps.g, alpha=alpha, Z1.params=Z1.params,Z2.params=Z2.params, 
                                mu=mu, beta.params=beta.params, eta.params = eta.params)
  
  # -- Outcome generation
  thr.weight <- NA
  thr.steps <- seq(0,1, step.size)
  betafn.true <- NA
  thresh_func <- switch(as.character(dg.spars), "weight-based"=cpp_weight_thresholding, "density-based"=cpp_density_thresholding)
  
  if(outcome %in% "prognostic") {
    if(dg.method %in% "single"){
      thr.weight=dg.thresh
      # Apply selected threshold to each ISN
      GE.thres <- data.frame(thresh_func(M=data.graph$GE, w=thr.weight, method = "trim"))
      # Compute graph features for each ISN
      GE.gvars <- data.frame(t(apply(GE.thres, 1, function(x) cpp_cc(x, p))))
      Xg <- unlist(GE.gvars)
    }else if(dg.method %in% "random"){
      thr.weight <- runif(n, min=dg.thresh[1], max=dg.thresh[2])
      GE.tmp <- lapply(1:nrow(data.graph$GE), function(x) thresh_func(matrix(data.graph$GE[x,], nrow=1), w=thr.weight[x], method = "trim"))
      GE.thres <- do.call(rbind,GE.tmp)
      GE.gvars <- data.frame(t(apply(GE.thres, 1, function(x) cpp_cc(x, p))))
      Xg <- unlist(GE.gvars)
    }else if(dg.method %in% "func"){
      betafn.true <- switch(dg.thresh, "flat"=rep(b1,length(thr.steps)),
                            "half-sine"=ifelse(thr.steps >0.5, 0, sin(thr.steps*pi*2)*b1*2), 
                            "sine"=ifelse(thr.steps >0.75, 0, sine_fn(thr.steps)))
      b1 <- switch(dg.thresh, "flat"=b1/5,"half-sine"=b1/5, "sine"=b1)
      GE.gvars.mat <- matrix(NA, ncol=length(thr.steps), nrow=n)
      for(t in 1:length(thr.steps)){
        GE.thres <- data.frame(thresh_func(data.graph$GE, thr.steps[t], method = "trim"))
        GE.gvars <- data.frame(t(apply(GE.thres, 1, function(x) cpp_cc(x, p=p))))
        GE.gvars.mat[,t] <- unlist(GE.gvars)
      }
      prod.beta.Xg <- GE.gvars.mat %*% diag(betafn.true)*step.size
      Xg <- rowSums(prod.beta.Xg)
    }
    
    X <- switch(as.character(setting), "uni"=0, "latent"=data.graph$Z[,2], "multi"=rnorm(n, mean=0, sd=0.25))
    b2 <- switch(as.character(setting), "uni"=0, "latent"=2, "multi"=2)
    xb <- b0 + Xg * b1 + X * b2
    Y <- rnorm(n, mean = xb, sd = eps.y)
    
    true.R2 = cor(Y, Xg)^2
    
  }else if(outcome %in% "diagnostic"){
      thr.weight = 0
      GE.thres <- data.frame(thresh_func(M=data.graph$GE, w=thr.weight, method = "trim"))
      # Compute graph features for each ISN
      GE.gvars <- data.frame(t(apply(GE.thres, 1, function(x) mean(WGCNA::clusterCoef(VecToSymMatrix(diag.entries=0,x, p))))))
      Xg <- unlist(GE.gvars)
    
      Y <- data.graph$Z[,1] + rnorm(n, mean = 0, sd = eps.y/10)
      X <- data.graph$Z[,2]
      
      true.R2 <- cor(Y, Xg)^2
  }else{
    stop("Parameter ouctome needs to be either prognostic or diagnostic!")
  }
  
  # --- Output
  out <- list("data"=data.frame("ID"=1:n, "Y"=Y, "Z"=data.graph$Z, "Xg"=Xg, "X"=X,
                                "GE"=data.graph$GE, "GE.noisy"=data.graph$GE.err,
                                "dg.method"=dg.method,"dg.threshold"=thr.weight, "true.R2"=true.R2),
              "fun"=data.frame("steps"=thr.steps,"betafn.true"=betafn.true),
              "network"=list("eta"=data.graph$eta, "eta.err"=data.graph$eta.err))
  return(out)   
}


# ====================== 02. Data analysis =====================================
analyse_data <- function(setting, df, network.model, n, p, b1, dg.thresh, dg.spars, k=5, step.size){
  #' Perform standard sparsification & flexible param approach
  #' setting=scn$setting; network.model=scn$network.model; df=data.iter; n=scn$n; p=scn$p; dg.thresh=scn$dg.thresh; k=5
  
  true.params <- data.frame("ID"= df$data$ID,
                           "DGMethod"=df$data$dg.method,
                           "Thresh"=df$data$dg.threshold,
                           "SparsMethod"="weight-based",
                           "ThreshMethod"="trim",
                           "Variable"="cc.uw")
  threshold.OPT <- data.frame(SparsMethod = c("weight-based", "density-based"), 
                          threshold.lo =c(0, .25),
                          threshold.up =c(.75, 1))
  threshold.AVG <- data.frame(SparsMethod = c("weight-based", "density-based"), 
                              threshold.lo =c(.1, .6),
                              threshold.up =c(.4, .9))
  dg.method <- ifelse(names(dg.thresh) %in% c("flat", "half-sine", "sine"), "func", names(dg.thresh) )
  adjust <- ifelse(setting %in% "uni", FALSE, TRUE)
  
  # Extract network data
  po = (p-1)*p/2                                                                  
  data.network <- df$data[,paste0("GE.noisy.",1:po)]
  df$data$fold <- cvFolds(length(unique(df$data$ID)), K=k)$which
  options(dplyr.summarise.inform = FALSE)
  
  # CC for threshold sequence
  data.gvars <- wrapperThresholding(df=data.network, msize=p, step.size=step.size)

  # Add outcome Y
  data.gvars <- merge(df$data[,c("ID","Y", "X", "fold")], data.gvars, by="ID") %>%
    mutate(Value=ifelse(is.nan(Value), NA, Value)) %>%
    mutate_at(vars(Thresh, Y, Value), as.numeric) %>%
    filter(Variable %in% "cc.uw") 
  
  # --  Oracle model
  # Data preparation for oracle model
  if(dg.method %in% c("single","random")){
    data.oracle <- data.gvars %>%
      filter(SparsMethod == true.params$SparsMethod & ThreshMethod == true.params$ThreshMethod &
               Variable == true.params$Variable) %>%
      group_by(Thresh) %>%
      mutate("true.t"=df$data$dg.threshold) %>%
      filter(Thresh == round(true.t,2)) %>%
      dplyr::select(!true.t)
  }else if(dg.method %in% "func"){
    mat.gvars <- data.gvars %>%
      filter(SparsMethod == true.params$SparsMethod & ThreshMethod == true.params$ThreshMethod &
               Variable == true.params$Variable) %>%
      dplyr::select(ID, Thresh, Value) %>%
      arrange(Thresh) %>%
      pivot_wider(names_from = "Thresh", values_from = "Value")
    
    prod.betaX.true <- rowSums(t(df$fun$betafn.true * t(mat.gvars[,-1])))
    
    data.oracle <- data.gvars %>%
      filter(SparsMethod == true.params$SparsMethod & ThreshMethod == true.params$ThreshMethod &
               Variable == true.params$Variable & Thresh == 0) %>%
      mutate(Value = prod.betaX.true*(b1/5))
  }
  
  data.oracle <- data.oracle %>%
    group_by(ThreshMethod, SparsMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) rbind(perform_AVG(dat=x, k=k, adjust=F, family = "gaussian"), 
                                              perform_AVG(dat=x, k=k, adjust=T, family = "gaussian")))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="Oracle",
           "Spline"=NA,
           Thresh=ifelse(length(dg.thresh)>1, "random", as.character(true.params$Thresh[1])))
  
  # --  Null model
  data.null <- df$data %>%
    dplyr::select(ID, fold, Y, X) %>%
    mutate(ThreshMethod = true.params$ThreshMethod,
           SparsMethod = true.params$SparsMethod,
           Variable = true.params$Variable,
           Value = mean(Y)) %>%
    group_by(ThreshMethod, SparsMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) rbind(perform_AVG(dat=x, k=k, adjust=F, family = "gaussian"),
                                              perform_AVG(dat=x, k=k, adjust=T, family = "gaussian")))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="Null",
           Thresh=as.character(Thresh)) %>% 
    slice(rep(1:n(), times = 2)) %>%
    mutate(SparsMethod = rep(c("weight-based", "density-based"), each=2))
    
  
  
  # -- Pick model with best RMSE
  data.OPT <- data.gvars  %>% 
    left_join(threshold.OPT, by = 'SparsMethod') %>%
    filter(Thresh >= threshold.lo & Thresh <= threshold.up) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) rbind(perform_OPT(dat=x, k=k, adjust=F, family = "gaussian"), 
                                              perform_OPT(dat=x, k=k, adjust=T, family = "gaussian")))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="OPT",
           Thresh=as.character(Thresh)) 
  
  
  # --  Average feature across threshold sequence
  data.AVG <- data.gvars %>%
    left_join(threshold.AVG, by = 'SparsMethod') %>%
    filter(Thresh >= threshold.lo & Thresh <= threshold.up) %>%
    group_by(SparsMethod, ThreshMethod, Variable, ID, Y, X, fold) %>%
    summarise("Value.avg"=mean(Value, na.rm=T)) %>%
    rename(Value=Value.avg) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) rbind(perform_AVG(dat=x, k=k, adjust=F, family = "gaussian"), 
                                              perform_AVG(dat=x, k=k, adjust=T, family = "gaussian")))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="AVG",
           Thresh=as.character(Thresh)) 
  
  
  # --  Functional data analysis approach
  data.FLEX <- data.gvars %>%
    arrange(Thresh) %>%
    mutate(Thresh=paste0("T_",Thresh),
           Y = as.numeric(as.character(Y))) %>%
    pivot_wider(values_from = Value, names_from = Thresh) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) rbind(perform_FLEX(data.fda=x, k=k, nodes=20, adjust=F, bs.type="ps", family = "gaussian"),
                                              perform_FLEX(data.fda=x, k=k, nodes=20, adjust=T, bs.type="ps", family = "gaussian")))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="FLEX",
           Thresh=as.character(Thresh))
  
  data.FLEX.coeff <- data.FLEX %>%
    dplyr::select(Variable, AnaMethod, Adjust, ThreshMethod, SparsMethod, Coef) %>%
    unnest(Coef) %>%
    reduce(data.frame) %>%
    `colnames<-`(c("Variable", "AnaMethod", "Adjust","ThreshMethod", "SparsMethod", 
                   "fda.thresh", "fda.est", "fda.se", "fda.sd.Xt.pred")) %>%
    merge(., df$fun, by.x="fda.thresh", by.y="steps")


  # -- Results
  out <- list()
  out$results <- data.frame(rbind(
    data.null[,c("AnaMethod", "Adjust", "SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.oracle[,c("AnaMethod","Adjust", "SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.OPT[,c("AnaMethod","Adjust", "SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.AVG[,c("AnaMethod","Adjust", "SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.FLEX[,c("AnaMethod","Adjust", "SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")]))
  out$results <- out$results %>%
    mutate(NetworkModel = network.model) %>%
    relocate(NetworkModel, .before=1)
  out$more$FLEX.coeff <- data.FLEX.coeff
  out$more$true.R2 <- df$data$true.R2[1]
  out$more$true.params <- true.params
  out$data <- data.gvars
  
  # Investigate reason for RMSE outliers 
  filename = paste0("data_",outcome, "_", setting, "_", network.model,"_n", n, "_p", p, "_b1", b1, "_DGT", dg.thresh, "_DGS", dg.spars)
  # if(any(data.FLEX[data.FLEX$SparsMethod %in% "weight-based",]$RMSE >8)){
  #   list_data <- list("pre_network"=df$network,"network"=data.network, "data"=data.gvars,"model"=data.FLEX, "coeffs"=data.FLEX.coeff)
  #   saveRDS(list_data, paste0(sim.path, filename , ".rds"))}
  
  return(out)
}


# =================== 03. Summarize & save scen results =============================
summarize_data <- function(results.sim, setting, outcome, n, p, q, network.model, mu, alpha0.params, alpha12.params, Z1.params, Z2.params, beta.params, eta.params, 
                           b0, b1, b2, eps.y, eps.g, epslevel.y, epslevel.g, dg.thresh, dg.spars, default.graph, filename, excel){
  #' Summarize results and save 
  # results.sim=results.sim; setting = scn$setting; outcome=scn$outcome; n=scn$n; p=scn$p; q=scn$q; network.model=scn$network.model;
  # alpha0.params=unlist(scn$alpha0.params, use.names = F); alpha12.params=unlist(scn$alpha12.params);
  # Z1.params=unlist(scn$Z1.params); Z2.params=unlist(scn$Z2.params); beta.params=unlist(scn$beta.params, use.names = F); eta.params=unlist(scn$eta.params);
  # b0=scn$b0; b1=scn$b1; b2=scn$b2; eps.y=scn$eps.y; eps.g=scn$eps.g; epslevel.y=scn$epslevel.y; epslevel.g=scn$epslevel.g
  # dg.thresh=scn$dg.thresh; default.graph = dnw.params$default.graph

  main.params <- data.frame(
    "iter" = iter,
    "n" = n,
    "q" = q,
    "p" = p,
    "setting" = setting,
    "outcome" = outcome,
    "network.model"=network.model,
    "dg.thresh" = names(dg.thresh),
    "dg.spars" = dg.spars,
    "beta.params" = paste0("beta(",beta.params[1],", ",beta.params[2],")"),
    "alpha0.params" = paste0("N(",alpha0.params[1],", ",alpha0.params[2]^2,")"),
    "alpha12.params" = paste0("U(",alpha12.params[1],", ",alpha12.params[2],")"),
    "Z1.params" = paste0("N(",Z1.params[1],", ",Z1.params[2]^2,")"),
    "Z2.params" = paste0("Ber(",Z2.params,")"),
    "b0" = b0,
    "b1" = b1,
    "b2" = b2,
    "epslevel.y" = epslevel.y,
    "epslevel.g" = epslevel.g,
    "eps.y" = eps.y,
    "eps.g" = eps.g
  )
  
  # -- Extraction of results 
  res <- list()
  
  # Overall performance comparison between methods 
  res$perf <- do.call(rbind, lapply(results.sim,  function(x) x[[1]])) %>%
    mutate(it=rep(1:iter, each=18)) %>%
    relocate(it, .before=1)
  
  # Coefficient function of the functional data approach
  res$more <- list()
  res$more$FLEX.coeff <- do.call(rbind, lapply(1:length(results.sim),  function(x) data.frame(iter=x, results.sim[[x]]$more$FLEX.coeff)))
  res$more$true.R2 <- do.call(rbind, lapply(1:length(results.sim),  function(x) data.frame(iter=x, results.sim[[x]]$more$true.R2)))
  true.R2 = mean(res$more$true.R2[,2], na.rm=T)    
  
  
  # -- Oracle 
  res.oracle <- res$perf %>%
    data.frame() %>%
    filter(AnaMethod %in% "Oracle") %>%
    group_by(AnaMethod, NetworkModel, Adjust, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.med"=median(RMSE, na.rm=T), 
              "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSE.est)  
  
  # -- Null 
  res.null <- res$perf %>%
    data.frame() %>%
    filter(AnaMethod %in% "Null") %>%
    group_by(AnaMethod, NetworkModel, Adjust, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.med"=median(RMSE, na.rm=T), 
              "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSE.est)  
  
  # -- Best RMSE 
  res.OPT <- res$perf %>%
    data.frame() %>%
    filter(AnaMethod %in% "OPT") %>%
    group_by(AnaMethod, NetworkModel, Adjust, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.med"=median(RMSE, na.rm=T), 
              "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSE.est)   
  
  res.OPT.freq <- res$perf %>%
    filter(AnaMethod %in% "OPT") %>%
    dplyr::count(Adjust, NetworkModel, SparsMethod, ThreshMethod, Thresh, Variable) %>%
    rename(count=n) %>%
    data.frame() %>%
    mutate(perc=round(count/iter*100,2)) %>%
    arrange(Variable, SparsMethod, ThreshMethod,desc(count)) 
  
  
  # -- Averaging over threshold sequence
  res.AVG <- res$perf %>%
    filter(AnaMethod %in% "AVG") %>%
    select(!Thresh) %>%
    group_by(AnaMethod, NetworkModel, Adjust, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.med"=median(RMSE, na.rm=T), 
              "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) %>%
    data.frame()
  
  
  # -- FLEX 
  res.FLEX <- res$perf %>%
    filter(AnaMethod %in% "FLEX") %>%
    select(!Thresh)  %>%
    group_by(AnaMethod, NetworkModel, Adjust, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.med"=median(RMSE, na.rm=T), 
              "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) %>%
    data.frame()
  
  res.FLEX.coeff <- res$more$FLEX.coeff %>%
    group_by(Adjust, SparsMethod, ThreshMethod, fda.thresh, Variable) %>%
    summarise("coeff.mean.ustd"=mean(fda.est, na.rm=T), 
              "coeff.median.ustd"=median(fda.est, na.rm=T), 
              "coeff.lo.ustd"=quantile(fda.est, 0.05, na.rm=T), 
              "coeff.up.ustd"=quantile(fda.est,0.95, na.rm=T),
              "coeff.mean.std"=mean(fda.est*fda.sd.Xt.pred, na.rm=T), 
              "coeff.median.std"=median(fda.est*fda.sd.Xt.pred, na.rm=T), 
              "coeff.lo.std"=quantile(fda.est*fda.sd.Xt.pred, 0.05, na.rm=T), 
              "coeff.up.std"=quantile(fda.est*fda.sd.Xt.pred,0.95, na.rm=T)) %>%
    data.frame()
  
  res.FLEX.func <- res$more$FLEX.coeff %>%
    group_by(Adjust, SparsMethod, ThreshMethod, fda.thresh, Variable) %>%
    summarise("aRMSE.ustd"=calc_rmse(fda.est, betafn.true),
              "aRMSE.std"=calc_rmse(fda.est*fda.sd.Xt.pred, betafn.true)) %>%
    data.frame()
  
  # -- Save results 
  main.params$true.R2 <- round_0(true.R2,4)
  main.params$default.graph.density <- round_0(edge_density(default.graph)*100,2)
  tbl_res <- data.frame(rbind(res.oracle, res.null, res.OPT, res.AVG, res.FLEX))
  list_results <- list("scenario"=main.params,
                       "results"=list("tbl_results"= cbind(main.params, tbl_res),
                                      "tbl_OPT_freq"=cbind(main.params, res.OPT.freq),
                                      "tbl_FLEX_coeff"=cbind(main.params, res.FLEX.coeff),
                                      "tbl_FLEX_func"=cbind(main.params, res.FLEX.func)),
                       "add"=list("graph_default"=default.graph),
                       "iters"=list("res"=cbind(main.params, res$perf),
                                    "fun"=cbind(main.params, res$more$FLEX.coeff)))
  
  
  saveRDS(list_results, here::here(sim.path, paste0(filename , ".rds")))  
  if(excel){write.xlsx(list_results, paste0(sim.path, filename, ".xlsx"), overwrite = T)}
  #return(list_results)
}



# ======================== 04. Concatenate simulation results ==================================
evaluate_scenarios <- function(sim.files, sim.path){
  res <- list()
  for(i in 1:length(sim.files)){
    list.tmp <- readRDS(here::here(sim.path, sim.files[i]))
    list.tmp$scenario$dg.thresh <- ifelse(str_detect(list.tmp$scenario$dg.thresh, "random"), "random", list.tmp$scenario$dg.thresh)
    
    g <- as.matrix(as_adjacency_matrix(list.tmp$add$graph_default))
    res_graph <- data.frame(cbind(list.tmp$scenario, g)) %>%
      group_by_at(colnames(list.tmp$scenario)) %>%
      nest()
    
    res[[i]] <- list("sim"=list.tmp$results$tbl_results, "fun"=list.tmp$results$tbl_FLEX_func, 
                     "tfreq"=list.tmp$results$tbl_OPT_freq, "funcoeff"=list.tmp$results$tbl_FLEX_coeff, 
                     "graph"=res_graph, "iters_res"=list.tmp$iters$res, "iters_fun"=list.tmp$iters$fun)
  }
  tbl_scens <- do.call(rbind, lapply(res, function(x) x[[1]]))
  tbl_funs <- do.call(rbind, lapply(res, function(x) x[[2]]))
  tbl_tfreq <- do.call(rbind, lapply(res, function(x) x[[3]]))
  tbl_funcoeff <- do.call(rbind, lapply(res, function(x) x[[4]]))
  tbl_graph <- do.call(rbind, lapply(res, function(x) x[[5]]))
  tbl_iters_res <- do.call(rbind, lapply(res, function(x) x[[6]]))
  tbl_iters_fun <- do.call(rbind, lapply(res, function(x) x[[7]]))
  
  write.csv(tbl_iters_res,paste0(here::here(sim.path, "/tbl_results_iters.csv")), row.names = FALSE)
  
  out <- list("tbl_scens"=tbl_scens, "tbl_funs"=tbl_funs, "tbl_tfreq"=tbl_tfreq, 
              "tbl_funcoeff"=tbl_funcoeff, "tbl_graph"=tbl_graph, 
              "tbl_iters_res"=tbl_iters_res, "tbl_iters_fun"=tbl_iters_fun)
  return(out)
}


# ======================== 05. Report simulation results ==================================
report_simresults <- function(sim.path, filename){
  
  filename=paste0("report_results_2022-10-11")
  rmarkdown::render(
    "src/report_main.Rmd",
    output_dir = sim.path,
    output_file = paste0(filename, ".html"))
  browseURL(file.path(paste0(sim.path, "/",filename, ".html")))
}
