# ============================================================================ #
# Author: MG
# Date: 19.04.2022
# Info: Main simulation functions 
# ============================================================================ #



# ================================== 00. Main  =====================================

simulate_pronet <- function(iter, n, p, q, b0, b1, dg.thresh, 
                            beta.params, alpha0.params, alpha12.params, Z1.params, Z2.params,
                            eps.y, eps.g, filename, excel){
  #' Given specific design parameters, performs a number of iterations and saves the result in a R object
  # iter=scn$iter; n=scn$n; p=scn$p; q=scn$q; b0=scn$b0; b1=scn$b1;
  # dg.thresh=scn$dg.thresh;
  # beta.params = scn$beta.params; alpha0.params = scn$alpha0.params; alpha12.params = scn$alpha12.params;
  # Z1.params = scn$Z1.params;  Z2.params = scn$Z2.params;
  # eps.y=scn$eps.y; eps.g=scn$eps.g

  # Preprocess
  beta.params = unlist(beta.params, use.names = F)
  alpha0.params = unlist(alpha0.params, use.names = F)
  alpha12.params = unlist(alpha12.params, use.names = F)
  Z1.params = unlist(Z1.params, use.names = F)
  Z2.params = unlist(Z2.params, use.names = F)

  # -- Setup default network
  dnw.params <- genDefaultNetwork(p, q, beta.params, alpha0.params, alpha12.params, Z1.params, Z2.params)
  
  # -- Data generation & analysis
  results.sim <- list()
  results.sim <- lapply(1:iter, function(x){
    data.iter <- generate_data(n = n, 
                               p = p, 
                               q = q,
                               alpha = dnw.params$alpha, 
                               mu = dnw.params$mu, 
                               eta.params = dnw.params$eta.params,
                               beta.params = beta.params,
                               Z1.params = Z1.params,
                               Z2.params = Z2.params,
                               b0 = b0,
                               b1 = b1,  
                               eps.y = eps.y, 
                               eps.g = eps.g, 
                               dg.thresh = dg.thresh)
    results.iter <- analyse_data(data.iter, 
                                 n = n, 
                                 p = p,
                                 dg.thresh = unlist(dg.thresh), 
                                 k = 5)
    return(results.iter)
  })

  # -- Summarize & save results
  summarize_data(results.sim, n=n, p=p, q=q, 
                 alpha0.params = alpha0.params, alpha12.params = alpha12.params, 
                 Z1.params=Z1.params, Z2.params=Z2.params,beta.params=beta.params, eta.params=eta.params, 
                 b0=b0, b1=b1, eps.y=eps.y, eps.g=eps.g, 
                 dg.thresh=unlist(dg.thresh), BA.graph = dnw.params$BA.graph, 
                 filename=filename, excel=excel)
}

# ============================ 01. Data generation =============================
generate_data <- function (n, p, q, mu, alpha, Z1.params, Z2.params, beta.params, eta.params, 
                           b0, b1, eps.y, eps.g, dg.thresh) {
  #' Data generation (code adapted and modified; initially from https://github.com/shanghongxie/Covariate-adjusted-network)
  # n=scn$n; p=scn$p; q=scn$q;
  # alpha=dnw.params$alpha; mu=dnw.params$mu; eta.params = dnw.params$eta.params;
  # beta.params = unlist(scn$beta.params); Z1.params = unlist(scn$Z1.params); Z2.params = unlist(scn$Z2.params);
  # b0=scn$b0; b1 = scn$b1;
  # eps.y=scn$eps.y; eps.g=scn$eps.g; dg.thresh=scn$dg.thresh

  # -- Individual-specific networks: Generate and analyse
  # Generate ISNs
  po = (p-1)*p/2    
  dg.method <- names(dg.thresh)
  dg.thresh <- unlist(dg.thresh)
  data.graph <- genIndivNetwork(n=n, p=p, q=q, alpha=alpha, Z1.params=Z1.params,Z2.params=Z2.params, 
                                mu=mu, beta.params=beta.params, eta.params = eta.params)
  GE <- abs(data.graph$GE)

  # Threshold ISN by a single cut-off for all indivs or select for each indiv a threshold within specified sequence
  step.size <- 0.02
  thr.weight <- NA
  thr.grid <- seq(0,1, step.size)
  betafn.true <- NA
  if(dg.method %in% "single"){
    thr.weight=dg.thresh
    # Apply selected threshold to each ISN
    GE.thres <- data.frame(cpp_weight_thresholding(M=GE, w=thr.weight, method = "trim"))
    # Compute graph features for each ISN
    GE.gvars <- data.frame(t(apply(GE.thres, 1, function(x) cpp_cc_func(x, p))))
    Xg <- unlist(GE.gvars)
  }else if(dg.method %in% "random"){
    thr.weight <- runif(n, min=dg.thresh[1], max=dg.thresh[2])
    GE.tmp <- lapply(1:nrow(GE), function(x) cpp_weight_thresholding(matrix(GE[x,], nrow=1), w=thr.weight[x], method = "trim"))
    GE.thres <- do.call(rbind,GE.tmp)
    GE.gvars <- data.frame(t(apply(GE.thres, 1, function(x) cpp_cc_func(x, p))))
    Xg <- unlist(GE.gvars)
  }else if(dg.method %in% "func"){
    if(dg.thresh %in% "half-sine"){
      betafn.true <- sin(thr.grid*pi)*b1*2
    }else if(dg.thresh %in% "sine"){
      betafn.true <- sin(thr.grid*pi*2)*b1*2
    }else if(dg.thresh %in% "flat"){
      betafn.true <- rep(b1*2,length(thr.grid))
    }
    
    GE.gvars.mat <- matrix(NA, ncol=length(thr.grid), nrow=n)
    for(t in 1:length(thr.grid)){
      GE.thres <- data.frame(cpp_weight_thresholding(GE, thr.grid[t], method = "trim"))
      GE.gvars <- data.frame(t(apply(GE.thres, 1, function(x) cpp_cc_func(x, p=p))))
      GE.gvars.mat[,t] <- unlist(GE.gvars)
    }
    Xg <- rowSums(GE.gvars.mat %*% betafn.true*step.size)
    b1 <- 1
  }
  
  # -- Outcome generation
  xb <- b0 + Xg * b1
  Y <- rnorm(n, mean = xb, sd = eps.y)
  
  true.R2 = cor(Y, Xg)^2
  
  # Generate noisy network data for data analyst
  GE.noisy.unscd <- abs(GE + sample(c(-1,1),size=n, replace=T) * rnorm(n, mean=0, sd=eps.g))
  GE.noisy <- t(apply(GE.noisy.unscd , 1, function(x) rowwise_scaling01(x)))
  out <- list("data"=data.frame("ID"=1:n,"Y"=Y, "Z"=data.graph$Z, "GE.fea"=Xg, 
                                "GE"=GE, "GE.noisy"=GE.noisy,
                   "dg.method"=dg.method,"dg.threshold"=thr.weight, "true.R2"=true.R2),
             "DG_funcform"=data.frame("grid"=thr.grid,"betafn.true"=betafn.true))
  return(out)   
}


# ====================== 02. Data analysis =====================================
analyse_data <- function(df, n, p, dg.thresh, k=5){
  #' Perform standard sparsification & flexible param approach
  #' df=data.iter; n=n; p=p; dg.thresh=scn$dg.thresh; k=5

  true.params = data.frame("ID"= df$data$ID,
                           "DGMethod"=df$data$dg.method,
                           "Thresh"=df$data$dg.threshold,
                           "SparsMethod"="weight-based",
                           "ThreshMethod"="trim",
                           "Variable"="cc.uw")
  threshold <- data.frame(SparsMethod = c("weight-based", "density-based"), 
                          threshold.lo =c(.1, .5),
                          threshold.up =c(.5, .9))
  
  # Extract network data
  po = (p-1)*p/2                                                                  
  data.network <- df$data[,paste0("GE.noisy.",1:po)]
  df$data$fold <- cvFolds(length(unique(df$data$ID)), K=k)$which
  options(dplyr.summarise.inform = FALSE)
  
  # CC for threshold sequence
  data.gvars <- wrapperThresholding(df=data.network, msize=p)

  # Add outcome Y
  data.gvars <- merge(df$data[,c("ID","Y", "fold")], data.gvars, by="ID") %>%
    mutate(Value=ifelse(is.nan(Value), NA, Value)) %>%
    mutate_at(vars(Thresh, Y, Value), as.numeric) %>%
    filter(Variable %in% "cc.uw") 
  
  # --  Oracle model
  data.oracle <- df$data %>%
    dplyr::select(ID, fold, Y, GE.fea) %>%
    rename(Value=GE.fea) %>%
    mutate(ThreshMethod = true.params$ThreshMethod,
           SparsMethod = true.params$SparsMethod,
           Variable = true.params$Variable) %>%
    group_by(ThreshMethod, SparsMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) evalLM(data.lm=x, k=k))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="Oracle",
           Thresh=ifelse(length(unique(true.params$Thresh))>1, "random", as.character(true.params$Thresh[1]))) 
  
  # --  Null model
  data.null <- df$data %>%
    dplyr::select(ID, fold, Y) %>%
    mutate(ThreshMethod = true.params$ThreshMethod,
           SparsMethod = true.params$SparsMethod,
           Variable = true.params$Variable,
           Value = mean(Y)) %>%
    group_by(ThreshMethod, SparsMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) evalLM(data.lm=x, k=k))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="Null",
           Thresh=ifelse(length(unique(true.params$Thresh))>1, "random", as.character(true.params$Thresh[1]))) 
  
  
  # -- Pick model with best RMSE
  data.bRMSE <- data.gvars  %>% 
    left_join(threshold, by = 'SparsMethod') %>%
    filter(Thresh >= threshold.lo & Thresh <= threshold.up) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) evalLM_dCV(data.lm=x,k=k))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="bRMSE",
           Thresh=as.character(Thresh)) 
  
  
  # --  Average feature across threshold sequence
  data.AVG <- data.gvars %>%
    left_join(threshold, by = 'SparsMethod') %>%
    filter(Thresh >= threshold.lo & Thresh <= threshold.up) %>%
    group_by(SparsMethod, ThreshMethod, Variable, ID, Y, fold) %>%
    summarise("Value.avg"=mean(Value, na.rm=T)) %>%
    rename(Value=Value.avg) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) evalLM(data.lm=x, k=k))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="AVG",
           Thresh=as.character(Thresh)) 
  
  
  # --  Functional data analysis approach
  data.FDA <- data.gvars %>%
    arrange(Thresh) %>%
    mutate(Thresh=paste0("T_",Thresh)) %>%
    pivot_wider(values_from = Value, names_from = Thresh) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(x) evalPFR(data.fda=x, k=k, bs.type="ps", nodes=20))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="FDA",
           Thresh=as.character(Thresh))
  
  data.FDA.coeff <- data.FDA %>%
    select(Variable, AnaMethod, ThreshMethod, SparsMethod, Coef) %>%
    unnest(Coef) %>%
    reduce(data.frame) %>%
    `colnames<-`(c("Variable", "AnaMethod", "ThreshMethod", "SparsMethod", 
                   "fda.thresh", "fda.est", "fda.se", "fda.sd.Xt.pred")) %>%
    merge(., df$DG_funcform, by.x="fda.thresh", by.y="grid")


  # -- Results
  out <- list()
  out$results <- data.frame(rbind(
    data.null[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.oracle[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.bRMSE[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.AVG[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.FDA[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")]))
  out$more$FDA.coeff <- data.FDA.coeff
  out$more$true.R2 <- df$data$true.R2[1]
  out$more$true.params <- true.params
  out$data <- data.gvars
  
  return(out)
}


# =================== 03. Summarize & save scen results =============================
summarize_data <- function(results.sim, n, p, q, mu, alpha0.params, alpha12.params, Z1.params, Z2.params, beta.params, eta.params, 
                           b0, b1, eps.y, eps.g, dg.thresh, BA.graph, filename, excel){
  #' Summarize results and save 
  # results.sim=results.sim; n=scn$n; p=scn$p; q=scn$q; alpha0.params=unlist(scn$alpha0.params, use.names = F); alpha12.params=unlist(scn$alpha12.params);
  # Z1.params=unlist(scn$Z1.params); Z2.params=unlist(scn$Z2.params); beta.params=unlist(scn$beta.params, use.names = F); eta.params=unlist(scn$eta.params);
  # b0=scn$b0; b1=scn$b1; eps.y=scn$eps.y; eps.g=scn$eps.g;
  # dg.thresh=scn$dg.thresh
  
  main.params <- list(
    iter = iter,
    n = n,
    q = q,
    p = p,
    dg.thresh = ifelse(names(dg.thresh) %in% "func", unlist(dg.thresh), names(dg.thresh)),
    beta.params = beta.params,
    alpha0.params = alpha0.params,
    alpha12.params = alpha12.params,
    Z1.params = Z1.params,
    Z2.params = Z2.params,
    b0 = b0,
    b1 = b1,
    eps.y = eps.y,
    eps.g = eps.g
  )
  
  # -- Extraction of results 
  res <- list()
  
  # Overall performance comparison between methods 
  res$perf <- do.call(rbind, lapply(results.sim,  function(x) x[[1]])) %>%
    mutate(iter=rep(1:iter, each=8))
  
  # Coefficient function of the functional data approach
  res$more <- list()
  res$more$FDA.coeff <- do.call(rbind, lapply(1:length(results.sim),  function(x) data.frame(iter=x, results.sim[[x]]$more$FDA.coeff)))
  res$more$true.R2 <- do.call(rbind, lapply(1:length(results.sim),  function(x) data.frame(iter=x, results.sim[[x]]$more$true.R2)))
  true.R2 = mean(res$more$true.R2[,2], na.rm=T)    
  
  
  # -- Oracle 
  res.oracle <- res$perf %>%
    data.frame() %>%
    filter(AnaMethod %in% "Oracle") %>%
    group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T),"R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T),"CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSE.est)  
  
  # -- Null 
  res.null <- res$perf %>%
    data.frame() %>%
    filter(AnaMethod %in% "Null") %>%
    group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T),"R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T),"CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSE.est)  
  
  # -- Best RMSE 
  res.bRMSE <- res$perf %>%
    data.frame() %>%
    filter(AnaMethod %in% "bRMSE") %>%
    group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T),"R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T),"CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSE.est)   
  
  res.bRMSE.freq <- res$perf %>%
    filter(AnaMethod %in% "bRMSE") %>%
    dplyr::count(SparsMethod, ThreshMethod, Thresh, Variable) %>%
    rename(count=n) %>%
    data.frame() %>%
    mutate(perc=round(count/iter*100,2)) %>%
    arrange(Variable, SparsMethod, ThreshMethod,desc(count)) 
  
  
  # -- Averaging over threshold sequence
  res.AVG <- res$perf %>%
    filter(AnaMethod %in% "AVG") %>%
    select(!Thresh) %>%
    group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T),"R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T),"CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) 
  
  
  # -- FDA 
  res.FDA <- res$perf %>%
    filter(AnaMethod %in% "FDA") %>%
    select(!Thresh)  %>%
    group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE, na.rm=T), "RMSE.lo"=quantile(RMSE, 0.05, na.rm=T), "RMSE.up"=quantile(RMSE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T),"R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T),"CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T)) 
  
  res.FDA.coeff <- res$more$FDA.coeff %>%
    group_by(SparsMethod, ThreshMethod, fda.thresh, Variable) %>%
    summarise("coeff.mean.ustd"=mean(fda.est, na.rm=T), 
              "coeff.lo.ustd"=quantile(fda.est, 0.05, na.rm=T), 
              "coeff.up.ustd"=quantile(fda.est,0.95, na.rm=T),
              "coeff.mean.std"=mean(fda.est*fda.sd.Xt.pred, na.rm=T), 
              "coeff.lo.std"=quantile(fda.est*fda.sd.Xt.pred, 0.05, na.rm=T), 
              "coeff.up.std"=quantile(fda.est*fda.sd.Xt.pred,0.95, na.rm=T)) %>%
    data.frame()
  
  res.FDA.func <- res$more$FDA.coeff %>%
    group_by(SparsMethod, ThreshMethod, fda.thresh, Variable) %>%
    summarise("aRMSE.ustd"=calc_rmse(fda.est, betafn.true),
              "aRMSE.std"=calc_rmse(fda.est*fda.sd.Xt.pred, betafn.true)) %>%
    data.frame()
  
  # -- Save results 
  main.params$true.R2 <- round_0(true.R2,4)
  main.params$BA.density <- round_0(edge_density(BA.graph)*100,2)
  main.params <- as_tibble(cbind(Setting = names(main.params), tibble(main.params))) %>%
    pivot_wider(names_from = Setting, values_from = main.params) %>%
    mutate(across(everything(), as.character)) %>%
    data.frame()
  tbl_res <- data.frame(rbind(res.oracle, res.null, res.bRMSE, res.AVG, res.FDA))
  list_results <- list("scenario"=main.params,
                       "tbl_results"= tbl_res,
                       "tbl_bRMSE_freq"=res.bRMSE.freq,
                       "tbl_FDA_coeff"=res.FDA.coeff,
                       "tbl_FDA_func"=res.FDA.func)
  
  
  saveRDS(list_results, paste0(sim.path, filename , ".rds"))  
  if(excel){write.xlsx(list_results, paste0(sim.path, filename, ".xlsx"), overwrite = T)}
  #return(list_results)
}


# ======================== 04. Report simulation results ==================================
report_simresults <- function(sim.path, filename){
  
  # Report results
  rmarkdown::render(
    "src/report_main.Rmd",
    params = list(output_dir=sim.path, scen.nr=scen.nr),
    output_dir = sim.path,
    output_file = paste0(filename, ".html"))
  browseURL(file.path(paste0(sim.path, filename, ".html")))
}


