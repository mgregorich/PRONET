# ============================================================================ #
# Author: MG
# Date: 19.04.2022
# Info: Main simulation functions 
# ============================================================================ #



# ================================== MAIN  =====================================

simulate_pronet <- function(main.params){
  #' Given specific design parameters, performs a number of iterations and saves the result in a R object
  
  # -- Setup default network
  dnw.params <- genDefaultNetwork(main.params, distr.params, BA.graph)
  
  # -- Data generation & analysis
  plan(multisession, workers=detectCores()/2)
  results.sim <- results.iter <-list()
  results.sim <- future_lapply(1:main.params$iter, function(x){
    data.iter <- generate_data(n=main.params$n, p=main.params$p, q=main.params$q,
                               alpha=dnw.params$alpha, mu=dnw.params$mu, eta.params = dnw.params$eta.params,
                               delta=main.params$delta, distr.params = distr.params, 
                               obeta0=main.params$obeta0, beta0=main.params$beta0,xbeta=main.params$xbeta,gbeta = main.params$gbeta,  
                               eps = list("y"=main.params$eps.y, "g"=main.params$eps.g), 
                               sthresh=main.params$sthresh)
    results.iter <- analyse_data(data.iter, true.params=true.params, tseq=da.thresh)
    return(results.iter)
  }, future.seed = 666)
  plan(sequential)  
  
  # -- Summarize & save results
  summarize_data(results.sim, distr.params, main.params)
}

# ============================ 01. Data generation =============================
generate_data <- function (n, p, q, mu, alpha, distr.params, eta.params, obeta0, 
                           delta, beta0, xbeta, gbeta, eps, sthresh) {
  # Data generation (code adapted and modified; initially from https://github.com/shanghongxie/Covariate-adjusted-network)
  # n=main.params$n; p=main.params$p; q=main.params$q; mu=main.params$mu; distr.params=distr.params;
  # alpha=main.params$alpha; delta=main.params$delta;obeta0=main.params$obeta0; beta0=main.params$beta0;xbeta=main.params$xbeta;
  # gbeta = main.params$gbeta; eps = list("y"=main.params$eps.y, "g"=main.params$eps.g);
  # sthresh=main.params$sthresh
  
  
  # -- Individual-specific network generation
  data.graph <- genIndivNetwork(n=n, p=p, q=q, alpha=alpha, 
                                distr.params=distr.params, eta.params = eta.params, delta, mu=mu)
  GE <- abs(data.graph$GE)
  GE.thres <- t(sapply(1:nrow(GE), function(x){
    if(length(sthresh)==1){
      thr.weight=sthresh
    }else{
      thr.weight=sample(sthresh,1, replace=T) 
    }
    mat <- VecToSymMatrix(0, side.entries = GE[x,], mat.size = p)
    res.mat <- data.frame(weightThresholding(mat, w=thr.weight, method = "trim")$adj)
    res <- res.mat[upper.tri(res.mat)]
  }))
  GE.fea <- data.frame(t(apply(GE.thres, 1, function(x) calcGraphFeatures(VecToSymMatrix(0, x, p)))))
  
  
  # -- Outcome generation
  Y=NULL
  for (i in 1:n) {
    xb=beta0+sum(GE.fea[i,1]*gbeta)
    Y[i]=rnorm(1,mean=xb,sd=eps$y)
  }
  GE.noisy = GE + matrix(rnorm(length(GE),0, sd=eps$g), nrow = 250)
  true.R2 = cor(Y, GE.fea[,1])^2
  
  
  df <- data.frame("Y"=Y, "X"=data.graph$X, "GE.fea"=GE.fea, "GE"=GE, "GE.noisy"=GE.noisy, "true.tau"=sthresh, "true.R2"=true.R2)
  return(df)   
}


# ====================== 02. Data analysis =====================================
analyse_data <- function(df, tseq, true.params, k=5){
  #' Apply standard sparsification & flexible param approach
  # df=data.test; tseq=da.thresh; k=5; 
  # true.params=list("SparsMethod"="weight-based", "ThreshMethod"="trim", "Thresh"=main.params$sthresh)
  
  # Extract network data
  data.network <- df[,paste0("GE.noisy.",1:po)]
  df$Subj <- 1:n
  df$fold <- cvFolds(length(unique(df$Subj)), K=k)$which
  options(dplyr.summarise.inform = FALSE)
  
  # CC for threshold sequence
  list.gvars <- lapply(1:nrow(data.network), function(x) data.frame("Subj"=x, wrapperThresholding(eweights=data.network[x,], msize=p, tseq=tseq)))
  data.gvars <- data.frame((do.call(rbind, list.gvars)))
  
  # Add outcome Y
  data.gvars <- merge(df[,c("Subj","Y", "fold")], data.gvars, by="Subj") %>%
    mutate(Value=ifelse(is.nan(Value), NA, Value)) %>%
    filter(Variable %in% "cc")
  
  # --  Oracle model
  data.oracle <- data.gvars  %>% 
    tibble() %>%
    filter(SparsMethod %in% true.params$SparsMethod & ThreshMethod %in% true.params$ThreshMethod & Thresh %in% true.params$Thresh) %>%
    select(!Thresh) %>%
    group_by(Variable, ThreshMethod, SparsMethod) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM(data.lm=df, k=k))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="Oracle",
           Thresh=true.params$Thresh)
  
  
  # -- Pick model with best RMSE
  data.bRMSE <- data.gvars  %>% 
    tibble() %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM_dCV(data.lm=df,k=k))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="bRMSE")
  
  
  # --  Average feature across threshold sequence
  data.AVG <- data.gvars %>%
    group_by(SparsMethod, ThreshMethod, Variable, Subj, Y, fold) %>%
    summarise("Value.avg"=mean(Value, na.rm=T)) %>%
    rename(Value=Value.avg) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM(data.lm=df, k=k))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="AVG")
  
  
  # --  Functional data analysis approach
  data.FDA <- data.gvars %>%
    pivot_wider(values_from = Value, names_from = Thresh) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalPFR(data.fda=df, k=k, tseq=tseq))) %>%
    unnest(res) %>%
    mutate("AnaMethod"="FDA")
  
  data.FDA.coeff <- data.FDA %>%
    select(Variable, AnaMethod, ThreshMethod, SparsMethod, Coef) %>%
    unnest(Coef) %>%
    reduce(data.frame) %>%
    `colnames<-`(c("Variable", "AnaMethod", "ThreshMethod", "SparsMethod", "fda.thresh", "fda.est", "fda.se"))
  
  # -- Results
  out <- list()
  out$results <- data.frame(rbind(
    data.oracle[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.bRMSE[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.AVG[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")],
    data.FDA[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE", "R2", "CS")]))
  out$more$FDA.coeff <- data.FDA.coeff
  out$more$true.R2 <- df$true.R2[1]
  out$data <- data.gvars
  
  return(out)
}


# =================== 03. Summarize & save results =============================
summarize_data <- function(results.sim, distr.params, main.params){
  
  # -- Extraction of results 
  res <- list()
  # Overall performance comparison between methods 
  res$perf <- do.call(rbind, lapply(results.sim,  function(x) x[[1]])) %>%
    mutate(iter=rep(1:iter, each=19))
  
  res$more <- list()
  # Coefficient function of the functional data approach
  res$more$FDA.coeff <- do.call(rbind, lapply(1:length(results.sim),  function(x) data.frame(iter=x, results.sim[[x]]$more$FDA.coeff)))
  res$more$true.R2 <- do.call(rbind, lapply(1:length(results.sim),  function(x) data.frame(iter=x, results.sim[[x]]$more$true.R2)))
  true.R2 = mean( res$more$true.R2[,2])
  
  
  # -- Oracle 
  res.oracle <- res$perf %>%
    data.frame() %>%
    filter(AnaMethod %in% "Oracle") %>%
    group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE), "RMSE.lo"=quantile(RMSE, 0.05), "RMSE.up"=quantile(RMSE, 0.95),
              "R2.est"=mean(R2),"R2.lo"=quantile(R2, 0.05), "R2.up"=quantile(R2, 0.95),
              "CS.est"=mean(CS),"CS.lo"=quantile(CS, 0.05), "CS.up"=quantile(CS, 0.95)) %>%
    data.frame() %>%
    arrange(RMSE.est)  
  
  
  # -- Best RMSE 
  res.bRMSE <- res$perf %>%
    data.frame() %>%
    filter(AnaMethod %in% "bRMSE") %>%
    group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE), "RMSE.lo"=quantile(RMSE, 0.05), "RMSE.up"=quantile(RMSE, 0.95),
              "R2.est"=mean(R2),"R2.lo"=quantile(R2, 0.05), "R2.up"=quantile(R2, 0.95),
              "CS.est"=mean(CS),"CS.lo"=quantile(CS, 0.05), "CS.up"=quantile(CS, 0.95)) %>%
    data.frame() %>%
    arrange(RMSE.est)   
  
  res.bRMSE.freq <- res$perf %>%
    filter(AnaMethod %in% "bRMSE" & SparsMethod %in% true.params$SparsMethod & ThreshMethod %in% true.params$ThreshMethod) %>%
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
    summarise("RMSE.est"=mean(RMSE), "RMSE.lo"=quantile(RMSE, 0.05), "RMSE.up"=quantile(RMSE, 0.95),
              "R2.est"=mean(R2),"R2.lo"=quantile(R2, 0.05), "R2.up"=quantile(R2, 0.95),
              "CS.est"=mean(CS),"CS.lo"=quantile(CS, 0.05), "CS.up"=quantile(CS, 0.95)) 
  
  
  # -- FDA 
  res.FDA <- res$perf %>%
    filter(AnaMethod %in% "FDA") %>%
    select(!Thresh)  %>%
    group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
    summarise("RMSE.est"=mean(RMSE), "RMSE.lo"=quantile(RMSE, 0.05), "RMSE.up"=quantile(RMSE, 0.95),
              "R2.est"=mean(R2),"R2.lo"=quantile(R2, 0.05), "R2.up"=quantile(R2, 0.95),
              "CS.est"=mean(CS),"CS.lo"=quantile(CS, 0.05), "CS.up"=quantile(CS, 0.95)) 
  
  res.FDA.coeff <- res$more$FDA.coeff %>%
    group_by(SparsMethod, ThreshMethod, fda.thresh, Variable) %>%
    summarise("coeff.mean"=mean(fda.est), "coeff.lo"=quantile(fda.est, 0.05), "coeff.up"=quantile(fda.est,0.95))
  
  
  # -- Save results 
  tbl_res <- data.frame(rbind(res.oracle,res.bRMSE,res.AVG, res.FDA))
  list_results <- list("main.params"=main.params,
                       "distr.params"=distr.params,
                       "tbl_results"= tbl_res,
                       "tbl_bRMSE_freq"=res.bRMSE.freq, 
                       "tbl_FDA_coeff"=res.FDA.coeff,
                       "data_example" =results.sim[[1]]$data,
                       "true.R2"=true.R2)
  filename <- paste0("sim_n",main.params$n,"_p",main.params$p,"_beta",distr.params$beta[1], "",
                     distr.params$beta[2],"_epsY",main.params$eps.y,".rds")
  saveRDS(list_results, paste0(sim.path, filename))  
  
  #return(list_results)
}


# ======================== 04. Report results ==================================
report_results <- function(sim.path){
  outfile <- paste0("report_n",main.params$n,"_p",main.params$p,"_beta",
                     distr.params$beta[1], "",distr.params$beta[2],"_epsY",main.params$eps.y,".html")
  # Report results
  rmarkdown::render(
    "main/report.Rmd",
    params = list(output_dir=sim.path),
    output_dir = sim.path,
    output_file = outfile)
  browseURL(file.path(paste0(sim.path,"Report_" ,Sys.Date(), ".html")))
  
}
