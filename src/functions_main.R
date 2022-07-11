# ============================================================================ #
# Author: MG
# Date: 19.04.2022
# Info: Main simulation functions 
# ============================================================================ #



# ================================== 00. Main  =====================================

simulate_pronet <- function(iter, n, p, q, b0, b1, dg.thresh, 
                            beta.params, alpha0.params, alpha12.params, X1.params, X2.params,
                            eps.y, eps.g, BA.graph, filename, excel){
  #' Given specific design parameters, performs a number of iterations and saves the result in a R object
  # iter=scn$iter; n=scn$n; p=scn$p; q=scn$q; b0=scn$b0; b1=scn$b1;
  # dg.thresh=scn$dg.thresh;
  # beta.params = scn$beta.params; alpha0.params = scn$alpha0.params; alpha12.params = scn$alpha12.params;
  # X1.params = scn$X1.params;  X2.params = scn$X2.params;
  # eps.y=scn$eps.y; eps.g=scn$eps.g; BA.graph=BA.graph

  # Preprocess
  beta.params = unlist(beta.params, use.names = F)
  alpha0.params = unlist(alpha0.params, use.names = F)
  alpha12.params = unlist(alpha12.params, use.names = F)
  X1.params = unlist(X1.params, use.names = F)
  X2.params = unlist(X2.params, use.names = F)
  if(is.list(dg.thresh)){dg.thresh = unlist(dg.thresh)}
  
  # -- Setup default network
  dnw.params <- genDefaultNetwork(p, q, BA.graph, beta.params, alpha0.params, alpha12.params, X1.params, X2.params)
  
  # -- Data generation & analysis
  results.sim <- list()
  plan(multisession, workers = detectCores()*.8)
  results.sim <- future_lapply(1:iter, function(x){
    data.iter <- generate_data(n = n, 
                               p = p, 
                               q = q,
                               alpha = dnw.params$alpha, 
                               mu = dnw.params$mu, 
                               eta.params = dnw.params$eta.params,
                               beta.params = beta.params,
                               X1.params = X1.params,
                               X2.params = X2.params,
                               b0 = b0,
                               b1 = b1,  
                               eps.y = eps.y, 
                               eps.g = eps.g, 
                               dg.thresh = dg.thresh)
    results.iter <- analyse_data(data.iter, 
                                 n = n, 
                                 p = p,
                                 dg.thresh = dg.thresh, 
                                 k = 5)
    return(results.iter)
  }, future.seed = 666)
  plan(sequential)    

  # -- Summarize & save results
  summarize_data(results.sim, n=n, p=p, q=q, 
                 alpha0.params = alpha0.params, alpha12.params = alpha12.params, 
                 X1.params=X1.params, X2.params=X2.params,beta.params=beta.params, eta.params=eta.params, 
                 b0=b0, b1=b1, eps.y=eps.y, eps.g=eps.g, 
                 dg.thresh=dg.thresh, BA.graph = BA.graph, 
                 filename=filename, excel=excel)
}

# ============================ 01. Data generation =============================
generate_data <- function (n, p, q, mu, alpha, X1.params, X2.params, beta.params, eta.params, 
                           b0, b1, eps.y, eps.g, dg.thresh) {
  #' Data generation (code adapted and modified; initially from https://github.com/shanghongxie/Covariate-adjusted-network)
  # n=scn$n; p=scn$p; q=scn$q;
  # alpha=dnw.params$alpha; mu=dnw.params$mu; eta.params = dnw.params$eta.params;
  # beta.params = unlist(scn$beta.params); X1.params = unlist(scn$X1.params); X2.params = unlist(scn$X2.params);
  # b0=scn$b0; b1 = scn$b1;
  # eps.y=scn$eps.y; eps.g=scn$eps.g; dg.thresh=unlist(scn$dg.thresh)

  # -- Individual-specific networks: Generate and analyse
  # Generate ISNs
  po = (p-1)*p/2                                                                  
  data.graph <- genIndivNetwork(n=n, p=p, q=q, alpha=alpha, X1.params=X1.params,X2.params=X2.params, mu=mu, beta.params=beta.params, eta.params = eta.params)
  GE <- abs(data.graph$GE)

  # Threshold ISN by a single cut-off for all indivs or select for each indiv a threshold within specified sequence
  if(length(dg.thresh)==1){
      thr.weight=rep(dg.thresh, n)
  }else{
      thr.weight=sample(dg.thresh, n, replace=T)
  }
  
  # Apply selected threshold to each ISN
  GE.thres <- matrix(NA, nrow=n, ncol = po)
  for(i in 1:n){
    mat <- VecToSymMatrix(0, side.entries = GE[i,], mat.size = p)
    res.mat <- data.frame(weightThresholding(mat, w=thr.weight[i], method = "trim")$adj)
    res <- res.mat[upper.tri(res.mat)]
    GE.thres[i,] <- res
  }
  # Compute graph features for each ISN
  GE.fea <- data.frame(t(apply(GE.thres, 1, function(x) calcGraphFeatures(VecToSymMatrix(0, x, p)))))
  
  
  # -- Outcome generation
  Y=NULL
  for (i in 1:n) {
    xb = b0 + sum(GE.fea[i,2] * b1)
    Y[i] = rnorm(1, mean=xb, sd=eps.y)
  }
  true.R2 = cor(Y, GE.fea[,1])^2
  GE.noisy = abs(scaling01(abs(GE + sample(c(-1,1),size=n, replace=T) * rnorm(n, mean=0, sd=eps.g))))
  
  df <- data.frame("Y"=Y, "X"=data.graph$X, "GE.fea"=GE.fea, "GE"=GE, "GE.noisy"=GE.noisy, 
                   "true.threshold"=thr.weight, "true.R2"=true.R2)
  return(df)   
}


# ====================== 02. Data analysis =====================================
analyse_data <- function(df, n, p, dg.thresh, k=5){
  #' Perform standard sparsification & flexible param approach
  #' df=data.iter; n=n; p=p; dg.thresh=unlist(scn$dg.thresh); k=5

  true.params = data.frame("Subj"= 1:nrow(df),
                           "Thresh"=df$true.threshold,
                           "SparsMethod"="weight-based",
                           "ThreshMethod"="trim",
                           "Variable"="cc")
  
  
  # Extract network data
  po = (p-1)*p/2                                                                  
  data.network <- df[,paste0("GE.noisy.",1:po)]
  df$Subj <- 1:n
  df$fold <- cvFolds(length(unique(df$Subj)), K=k)$which
  options(dplyr.summarise.inform = FALSE)
  
  # CC for threshold sequence
  list.gvars <- lapply(1:nrow(data.network), 
                       function(x) data.frame("Subj"=x, wrapperThresholding(eweights=data.network[x,], 
                                                                            msize=p, 
                                                                            tseq=seq(0,1,0.02))))
  data.gvars <- data.frame(do.call(rbind, list.gvars, quote=TRUE))
  
  
  # Add outcome Y
  data.gvars <- merge(df[,c("Subj","Y", "fold")], data.gvars, by="Subj") %>%
    mutate(Value=ifelse(is.nan(Value), NA, Value)) %>%
    filter(Variable %in% "cc")
  
  # --  Oracle model
  data.oracle <- df %>%
    dplyr::select(Subj, fold, Y, GE.fea.cc) %>%
    rename(Value=GE.fea.cc) %>%
    mutate(ThreshMethod = true.params$ThreshMethod,
           SparsMethod = true.params$SparsMethod,
           Variable = true.params$Variable) %>%
    group_by(ThreshMethod, SparsMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM(data.lm=df, k=k))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="Oracle",
           Thresh=ifelse(length(unique(true.params$Thresh))>1, "random", as.character(true.params$Thresh[1]))) 
  
  # --  Null model
  data.null <- df %>%
    dplyr::select(Subj, fold, Y) %>%
    mutate(ThreshMethod = true.params$ThreshMethod,
           SparsMethod = true.params$SparsMethod,
           Variable = true.params$Variable,
           Value = mean(Y)) %>%
    group_by(ThreshMethod, SparsMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM(data.lm=df, k=k))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="Null",
           Thresh=ifelse(length(unique(true.params$Thresh))>1, "random", as.character(true.params$Thresh[1]))) 
  
  
  # -- Pick model with best RMSE
  data.bRMSE <- data.gvars  %>% 
    filter(Thresh >=0.1 & Thresh <= 0.5) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM_dCV(data.lm=df,k=k))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="bRMSE",
           Thresh=as.character(Thresh)) 
  
  
  # --  Average feature across threshold sequence
  data.AVG <- data.gvars %>%
    filter(Thresh >=0.1 & Thresh <= 0.5) %>%
    group_by(SparsMethod, ThreshMethod, Variable, Subj, Y, fold) %>%
    summarise("Value.avg"=mean(Value, na.rm=T)) %>%
    rename(Value=Value.avg) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM(data.lm=df, k=k))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="AVG",
           Thresh=as.character(Thresh)) 
  
  
  # --  Functional data analysis approach
  data.FDA <- data.gvars %>%
    filter(Thresh >=0.1 & Thresh <= 0.9) %>%
    pivot_wider(values_from = Value, names_from = Thresh) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalPFR(data.fda=df, k=k))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="FDA",
           Thresh=as.character(Thresh))
  
  data.FDA.coeff <- data.FDA %>%
    select(Variable, AnaMethod, ThreshMethod, SparsMethod, Coef) %>%
    unnest(Coef) %>%
    reduce(data.frame) %>%
    `colnames<-`(c("Variable", "AnaMethod", "ThreshMethod", "SparsMethod", "fda.thresh", "fda.est", "fda.se")) 

  
  # -- Results
  out <- list()
  out$results <- data.frame(rbind(
    data.null[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.oracle[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.bRMSE[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.AVG[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.FDA[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")]))
  out$more$FDA.coeff <- data.FDA.coeff
  out$more$true.R2 <- df$true.R2[1]
  out$more$true.params <- true.params
  out$data <- data.gvars
  
  return(out)
}


# =================== 03. Summarize & save scen results =============================
summarize_data <- function(results.sim, n, p, q, mu, alpha0.params, alpha12.params, X1.params, X2.params, beta.params, eta.params, 
                           b0, b1, eps.y, eps.g, dg.thresh, BA.graph, filename, excel){
  #' Summarize results and save 
  # results.sim=results.sim; n=scn$n; p=scn$p; q=scn$q; alpha0.params=unlist(scn$alpha0.params, use.names = F); alpha12.params=unlist(scn$alpha12.params);
  # X1.params=unlist(scn$X1.params); X2.params=unlist(scn$X2.params); beta.params=unlist(scn$beta.params, use.names = F); eta.params=unlist(scn$eta.params);
  # b0=scn$b0; b1=scn$b1; eps.y=scn$eps.y; eps.g=scn$eps.g;
  # dg.thresh=scn$dg.thresh
  
  main.params <- list(
    iter = iter,
    n = n,
    q = q,
    p = p,
    dg.thresh = dg.thresh,
    beta.params = beta.params,
    alpha0.params = alpha0.params,
    alpha12.params = alpha12.params,
    X1.params = X1.params,
    X2.params = X2.params,
    b0 = b0,
    b1 = b1,
    eps.y = eps.y,
    eps.g = eps.g
  )
  
  # -- Extraction of results 
  res <- list()
  
  # Overall performance comparison between methods 
  res$perf <- do.call(rbind, lapply(results.sim,  function(x) x[[1]])) %>%
    mutate(iter=rep(1:iter, each=20))
  
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
    count(SparsMethod, ThreshMethod, Thresh, Variable) %>%
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
    summarise("coeff.mean"=mean(fda.est, na.rm=T), "coeff.lo"=quantile(fda.est, 0.05, na.rm=T), "coeff.up"=quantile(fda.est,0.95, na.rm=T)) %>%
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
                       "tbl_FDA_coeff"=res.FDA.coeff)
  
  
  saveRDS(list_results, paste0(sim.path, filename , ".rds"))  
  if(excel){write.xlsx(list_results, paste0(sim.path, filename, ".xlsx"), overwrite = T)}
  #return(list_results)
}


# ======================== 04. Report scen results ==================================
report_scenarios <- function(sim.path, scen.nr, filename){

  # Report results
  rmarkdown::render(
    "src/report_aux.Rmd",
    params = list(output_dir=sim.path, scen.nr=scen.nr),
    output_dir = sim.path,
    output_file = paste0(filename, ".html"))
  browseURL(file.path(paste0(sim.path, filename, ".html")))
}

# ======================== 05. Report simulation results ==================================
report_simresults <- function(sim.path, filename){
  
  # Report results
  rmarkdown::render(
    "src/report_main.Rmd",
    params = list(output_dir=sim.path, scen.nr=scen.nr),
    output_dir = sim.path,
    output_file = paste0(filename, ".html"))
  browseURL(file.path(paste0(sim.path, filename, ".html")))
}


