# ============================================================================ #
# Author: MG
# Date: 19.04.2022
# Info: Main simulation functions 
# ============================================================================ #



# ================================== 00. Simulate scenario  =====================================

simulate_scenario <- function(scn){
  #' Run each scenario with R simulation replicates

  sourceCpp(here::here("src","utils.cpp"))
  
  data_gen_thresh_short <- ifelse(scn$data_gen_thresh %in% "weight-based", "w", "d")
  filename <- paste0("sim_i", scn$iter, "_",scn$setting, "_",scn$data_gen_feature,"_n",scn$n,"_p",scn$p, 
                     "_DG",data_gen_thresh_short,"-",scn$data_gen_mech, 
                     "_1b",scn$b1,"_2b",scn$b2,"_epsY",scn$epslevel_y,"_epsG",scn$epslevel_g)    
  
  # -- Data generation & analysis
  results.sim <- list()
  results.sim <- lapply(1:scn$iter, function(x){
    data.iter <- generate_data(setting = scn$setting,
                               n = scn$n, 
                               data_gen_feature = scn$data_gen_feature,
                               p = scn$p, 
                               b0 = scn$b0,
                               b1 = scn$b1,  
                               b2 = scn$b2,
                               eps_y = scn$eps_y, 
                               eps_g = scn$eps_g, 
                               data_gen_thresh = scn$data_gen_thresh,
                               data_gen_mech = scn$data_gen_mech)
    results.iter <- analyse_data(setting = scn$setting,
                                 df = data.iter, 
                                 n = scn$n, 
                                 p = scn$p,
                                 b1 = scn$b1,
                                 data_gen_feature = scn$data_gen_feature,
                                 data_gen_thresh = scn$data_gen_thresh,
                                 data_gen_mech = scn$data_gen_mech)
    return(results.iter)
  })
  
  
  # -- Summarize & save results
  summarize_data(results.sim, setting=scn$setting, data_gen_feature=scn$data_gen_feature, n=scn$n, p=scn$p,
                 b0=scn$b0, b1=scn$b1, b2=scn$b2, eps_y=scn$eps_y, eps_g=scn$eps_g, epslevel_y=scn$epslevel_y, epslevel_g=scn$epslevel_g, 
                 data_gen_thresh=scn$data_gen_thresh, data_gen_mech=scn$data_gen_mech, filename=filename)
}

# ============================ 01. Data generation =============================
generate_data <- function (setting, n, data_gen_feature, p, b0, b1, b2, eps_y, eps_g, data_gen_thresh, data_gen_mech) {
  #' Data generation using a plasmode simulation design
  #' Data from the ABIDE preprocessed initiative is used: http://preprocessed-connectomes-project.org/abide/
  #' Outcome-generating mechanism based on underlying weighting function

  # -- Individual-specific networks generation
  # Generate ISNs
  po = (p-1)*p/2 
  data_gen_mech = as.character(data_gen_mech)
  
  # Subsample selection
  data_abide <- readRDS(here::here("data", "data_abide_plasmode.rds"))
  set_ids <- unique(data_abide$clin$ID)
  subset_ids <- sample(set_ids, n)
  
  data.clin <- data_abide$clin[data_abide$clin$ID %in% subset_ids,]
  data.gvars <- data_abide$gvars[data_abide$gvars$ID %in% subset_ids,] %>%
    arrange(ID, data_ana_thresh, data_ana_t)
  data.graph <- data_abide$graph[data_abide$graph$ID %in% subset_ids,] %>%
    dplyr::select(!ID) %>%
    as.matrix() 
  colnames(data.graph) <- NULL
  data.graph <- apply(data.graph,2,as.numeric)
  
  # -- Outcome generation
  data_gen_t <- NA
  grid_fun <- seq(0,1, step_size)
  betafun_true <- NA
  
  if(data_gen_mech %in% "universal"){
    data_gen_t <- ifelse(data_gen_thresh %in% "weight-based", data_gen_mech_universal_par, 1-2*data_gen_mech_universal_par)
    Xg <- data.gvars %>% filter(data_ana_thresh %in% data_gen_thresh & data_ana_t==data_gen_t & variable %in% data_gen_feature) %>% pull(value)
  }else if(data_gen_mech %in% "random"){
    if(data_gen_thresh %in% "weight-based"){
      data_gen_t <- data.frame("ID"=data.clin$ID, "true.thresh"=round(runif(n, min=data_gen_mech_random_pars[1], max=data_gen_mech_random_pars[2]),2))
    }else{
      data_gen_t <- data.frame("ID"=data.clin$ID, "true.thresh"=round(runif(n, min=1-data_gen_mech_random_pars[2], max=1-data_gen_mech_random_pars[1]),2))
      }
                  
    Xg <- full_join(data.gvars, data_gen_t, by="ID") %>%
      filter(data_ana_thresh %in% data_gen_thresh & variable %in% data_gen_feature) %>%
      mutate(dis=abs(data_ana_t-true.thresh)) %>%
      group_by(ID) %>%
      arrange(dis) %>%
      slice(1) %>%
      pull(value)
    data_gen_t <- data_gen_t$true.thresh
  }else if(data_gen_mech %in% c("flat","early peak","arc")){
    betafun_true <- switch(data_gen_mech, 
                           "flat"=rep(2,length(grid_fun)),
                           "early peak"=ifelse(grid_fun >0.5, 0, sin(grid_fun*pi*2)*3), 
                           "arc"=sin(grid_fun*pi)*3)
    b1 <- switch(data_gen_mech, "flat"=b1,"early peak"=b1, "arc"=b1)
    GE.gvars.mat <- data.gvars %>% 
      filter(data_ana_thresh %in% data_gen_thresh & variable %in% data_gen_feature) %>% 
      mutate(data_ana_t=paste0("T_", data_ana_t)) %>% 
      pivot_wider(names_from = data_ana_t, values_from = value) %>% 
      dplyr::select(paste0("T_", seq(0,1,step_size)))
    prod.beta.Xg <- t(t(GE.gvars.mat) * betafun_true)*step_size
    Xg <- rowSums(prod.beta.Xg)
  }

  if(setting=="uni"){
    b2 <- 0
    X <- 0
    xb <- b0 + Xg * b1 
    Y <- rnorm(n, mean = xb, sd = eps_y)  
  }else if(setting=="multi"){
    b2 <- 15
    X <- cbind(scaling01(data.clin$age), to_numeric(data.clin$diagnosis))
    xb <- b0 + Xg * b1 + X %*% c(b2, b2/3)
    Y <- rnorm(n, mean = xb, sd = eps_y)        
  }
  
  
  ## Noise settings
  data.graph.noisy <- matrix(0, nrow = n, ncol = po)
  
  if(eps_g!=0){
    scaling_factor <- round(runif(n, min=0, max=1),2)
    for(i in 1:n){
      GE <- VecToSymMatrix(diag.entries = 1, side.entries= data.graph[i,], mat.size = p)
      GE.noisy <- noisecor(cormat = GE, epsilon = eps_g, uidim = 2, delta=scaling_factor[i])
      data.graph.noisy[i,] <- GE.noisy[upper.tri(GE.noisy)]
    }    
  }else{
    scaling_factor <- 0
    data.graph.noisy <- data.graph
  }

  # Only consider positive connections
  data.graph.noisy[data.graph.noisy<0] <- 0
  data.graph.noisy[data.graph.noisy>1] <- 1
  
  # --- Output
  out <- list("data"=data.frame("ID"=data.clin$ID, "Y"=Y, "Xg"=Xg, "X"=X,
                                "GE"=data.graph, "GE.noisy"=data.graph.noisy,
                                "data_gen_thresh"=data_gen_thresh, "data_gen_mech"=data_gen_mech, "data_gen_t"=data_gen_t,
                                "noise_scaling"=scaling_factor),
              "fun"=data.frame("grid"=grid_fun,"betafun_true"=betafun_true, "b1"=b1))
  return(out)   
}


# ====================== 02. Data analysis =====================================
analyse_data <- function(df, setting, n, p, b1, data_gen_feature, data_gen_thresh, data_gen_mech, k=5){
  #' Perform ad-hoc competitors & flexible param approach

  true.params <- data.frame("ID"= df$data$ID,
                            "data_gen_thresh"=df$data$data_gen_thresh,
                            "data_gen_mech"=df$data$data_gen_mech,
                            "data_gen_t"=df$data$data_gen_t,
                            "variable"=data_gen_feature)
  # Extract network data
  po <- (p-1)*p/2 
  adjusted_ana <- ifelse(setting %in% "uni", FALSE, TRUE)
  data_gen_mech <- as.character(data_gen_mech)
  data.network <- df$data[,paste0("GE.noisy.",1:po)]
  df$data$fold <- sample(rep(1:5, ceil(n/k)))[1:n]
  options(dplyr.summarise.inform = FALSE)
  
  # CC for threshold sequence
  data.res <- wrapper_thresholding(dat=data.network, set_ids=unique(df$data$ID), 
                                   mdim=p, step_size = step_size, graph_feature=data_gen_feature)
  
  # Add outcome Y
  if(setting=="uni"){
    data.gvars <- merge(df$data[,c("ID","Y", "X", "fold")], data.res, by="ID") %>%
      mutate(value=ifelse(is.nan(value), NA, value)) %>%
      mutate_at(vars(data_ana_t, Y, value), as.numeric) %>%
      filter(variable %in% data_gen_feature) 
  }else{
    data.gvars <- merge(df$data[,c("ID","Y", "X.1", "X.2","fold")], data.res, by="ID") %>%
      mutate(value=ifelse(is.nan(value), NA, value)) %>%
      mutate_at(vars(data_ana_t, Y, value), as.numeric) %>%
      filter(variable %in% data_gen_feature)    
  }

  
  # --  Oracle model
  # Data preparation for oracle model
  if(data_gen_mech %in% c("universal","random")){
    data.oracle <- data.gvars %>%
      filter(data_ana_thresh == true.params$data_gen_thresh &
               variable == true.params$variable) %>%
      group_by(data_ana_t) %>%
      mutate("true.t"=df$data$data_gen_t) %>%
      filter(data_ana_t == round(true.t,2)) %>%
      dplyr::select(!true.t) 
  }else if(data_gen_mech %in% c("flat", "early peak", "arc")){
    mat.gvars <- data.gvars %>%
      filter(data_ana_thresh == true.params$data_gen_thresh &
               variable == true.params$variable) %>%
      dplyr::select(ID, data_ana_t, value) %>%
      arrange(data_ana_t) %>%
      pivot_wider(names_from = "data_ana_t", values_from = "value") %>%
      dplyr::select(!ID)
    prod.beta.Xg <- t(t(mat.gvars) * df$fun$betafun_true)*step_size
    Xg <- rowSums(prod.beta.Xg)
    
    data.oracle <- data.gvars %>%
      filter(data_ana_thresh == true.params$data_gen_thresh &
               variable == true.params$variable & data_ana_t == 0) %>%
      mutate(value = Xg)
  }
  
  unused_data_gen_thresh <- ifelse(data_gen_thresh %in% "density-based", "weight-based", "density-based")
  data.oracle <- data.oracle %>%
    group_by(data_ana_thresh, variable) %>%
    nest() %>%
    mutate(res.oracle=lapply(data, function(x) perform_AVG(dat=x, k=k, adjust=adjusted_ana, family = "gaussian"))) %>%
    unnest(res.oracle) %>%
    mutate("data_ana_model"="Oracle") %>%
    bind_rows(. , expand.grid("data_ana_model"="Oracle", "adjust"=adjusted_ana, "data_ana_thresh"=unused_data_gen_thresh, 
                              "data_ana_t"=NA, "variable"=data_gen_feature,"RMSPE"=NA, "R2"=NA, "CS"=NA))
  
  
  # --  Null model
  data.null <- df$data %>%
   # dplyr::select(ID, fold, Y, X) %>%
    mutate(data_ana_thresh = data_gen_thresh,
           variable = data_gen_feature,
           value = mean(Y)) %>%
    group_by(data_ana_thresh, variable) %>%
    nest() %>%
    mutate(res.null=lapply(data, function(x) perform_AVG(dat=x, k=k, adjust=adjusted_ana, family = "gaussian"))) %>%
    unnest(res.null) %>%
    mutate("data_ana_model"="Null",
           "data_ana_t"=NA) %>% 
    slice(rep(1:n(), times = 2)) %>%
    mutate(data_ana_thresh = rep(c("weight-based", "density-based"), each=1))
  
  
  
  # -- OPT model according to RMSPE
  threshold_OPT <- data.frame(data_ana_thresh = c("weight-based", "density-based"), threshold_lo =c(0, .5), threshold_up =c(.5, 1))
  data.OPT <- data.gvars  %>% 
    left_join(threshold_OPT, by = 'data_ana_thresh') %>%
    filter(data_ana_t >= threshold_lo & data_ana_t <= threshold_up) %>%
    group_by(data_ana_thresh, variable) %>%
    nest() %>%
    mutate(res.OPT=lapply(data, function(x) perform_OPT(dat=x, k=k, adjust=adjusted_ana, family = "gaussian"))) %>%
    unnest(res.OPT, keep_empty = T) %>%
    mutate("data_ana_model"="OPT",
           "data_ana_t"=as.character(data_ana_t)) 
  
  
  # --  AVG model: average feature across threshold sequence
  threshold_AVG <- data.frame(data_ana_thresh = c("weight-based", "density-based"),  threshold_lo =c(.1, .6), threshold_up =c(.4, .9))
  data.AVG <- data.gvars %>%
    left_join(threshold_AVG, by = 'data_ana_thresh') %>%
    filter(data_ana_t >= threshold_lo & data_ana_t <= threshold_up) %>%
    group_by(data_ana_thresh, variable, ID, Y, X, fold) %>%
    summarise("value.avg"=mean(value, na.rm=T)) %>%
    dplyr::rename(value=value.avg) %>%
    group_by(data_ana_thresh, variable) %>%
    nest() %>%
    mutate(res.AVG=lapply(data, function(x) perform_AVG(dat=x, k=k, adjust=adjusted_ana, family = "gaussian"))) %>%
    unnest(res.AVG, keep_empty = T) %>%
    mutate("data_ana_model"="AVG",
           "data_ana_t"=as.character(data_ana_t)) 
  
  
  # --  FLEX model: flexible parametrization of graph-theoretical feature sequence
  data.FLEX <- data.gvars %>%
    arrange(data_ana_t) %>%
    mutate(data_ana_t=paste0("T_",data_ana_t),
           Y = as.numeric(as.character(Y)),
           SM=data_ana_thresh) %>%
    pivot_wider(values_from = value, names_from = data_ana_t) %>%
    group_by(data_ana_thresh, variable) %>%
    nest() %>%
    mutate(res.FLEX=lapply(data, function(x) perform_FLEX(data.fda=x, k=k, adjust=adjusted_ana, family = "gaussian"))) %>%
    unnest(res.FLEX) %>%
    mutate("data_ana_model"="FLEX",
           "data_ana_t"=as.character(data_ana_t))
  
  data.FLEX.coeff <- data.FLEX %>%
    dplyr::select(variable, data_ana_model, adjust, data_ana_thresh, coef) %>%
    unnest(coef) %>%
    reduce(data.frame) %>%
    `colnames<-`(c("variable", "data_ana_model", "adjust", "data_ana_thresh", 
                   "fda.thresh", "fda.est", "fda.se", "fda.sd.Xt.pred")) %>%
    merge(., df$fun, by.x="fda.thresh", by.y="grid")
  
  
  # -- Results
  res <- data.frame(rbind(
    data.null[,c("data_ana_model", "adjust", "data_ana_thresh","data_ana_t","variable","RMSPE", "R2", "CS")],
    data.oracle[,c("data_ana_model","adjust", "data_ana_thresh", "data_ana_t","variable","RMSPE", "R2", "CS")],
    data.OPT[,c("data_ana_model","adjust", "data_ana_thresh", "data_ana_t","variable","RMSPE", "R2", "CS")],
    data.AVG[,c("data_ana_model","adjust", "data_ana_thresh", "data_ana_t","variable","RMSPE", "R2", "CS")],
    data.FLEX[,c("data_ana_model","adjust", "data_ana_thresh", "data_ana_t","variable","RMSPE", "R2", "CS")])) %>%
    group_by(data_ana_thresh, adjust) %>%
    mutate(relRMSPE = RMSPE/min(RMSPE, na.rm=T))
  out <- list()
  out$results <- res
  out$more$FLEX_coeff <- data.FLEX.coeff
  out$more$true.params <- true.params
  out$data <- data.gvars
  out$network <- df$network
  
  return(out)
}


# =================== 03. Summarize & save scen results =============================
summarize_data <- function(results.sim, setting,  data_gen_feature, n, p,   
                           b0, b1, b2, eps_y, eps_g, epslevel_y, epslevel_g, data_gen_thresh, data_gen_mech, filename){
  #' Summarize results across simulation replicates and save 
  
  main.params <- data.frame(
    "iter" = iter,
    "n" = n,
    "p" = p,
    "setting" = setting,
    "data_gen_feature" = data_gen_feature,
    "data_gen_thresh" = data_gen_thresh,
    "data_gen_mech" = data_gen_mech,
    "epslevel_y" = epslevel_y,
    "epslevel_g" = epslevel_g,
    "eps_y" = eps_y,
    "eps_g" = eps_g,
    "b0" = b0,
    "b1" = b1,
    "b2" = b2
  )
  
  # -- Extraction of results 
  res <- list()
  
  # Overall performance comparison between methods 
  size_tbl <- nrow(results.sim[[1]]$results)
  res$perf <- do.call(rbind, lapply(results.sim,  function(x) x[[1]])) %>%
    data.frame() %>%
    mutate(it=rep(1:iter, each=size_tbl)) %>%
    relocate(it, .before=1)
  
  # Coefficient function of the functional data approach
  res$more <- list()
  res$more$FLEX_coeff <- do.call(rbind, lapply(1:length(results.sim),  function(x) data.frame(iter=x, results.sim[[x]]$more$FLEX_coeff)))
  
  # -- Oracle 
  res_oracle <- res$perf %>%
    data.frame() %>%
    filter(data_ana_model %in% "Oracle") %>%
    group_by(data_ana_model, adjust, data_ana_thresh, variable) %>%
    summarise("RMSPE.est"=mean(RMSPE, na.rm=T), "RMSPE.med"=median(RMSPE, na.rm=T), 
              "RMSPE.lo"=quantile(RMSPE, 0.05, na.rm=T), "RMSPE.up"=quantile(RMSPE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T),
              "relRMSPE.est"=mean(relRMSPE, na.rm=T), "relRMSPE.med"=median(relRMSPE, na.rm=T), 
              "relRMSPE.lo"=quantile(relRMSPE, 0.05, na.rm=T), "relRMSPE.up"=quantile(relRMSPE, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSPE.est)  
  
  # -- Null 
  res_null <- res$perf %>%
    data.frame() %>%
    filter(data_ana_model %in% "Null") %>%
    group_by(data_ana_model, adjust, data_ana_thresh, variable) %>%
    summarise("RMSPE.est"=mean(RMSPE, na.rm=T), "RMSPE.med"=median(RMSPE, na.rm=T), 
              "RMSPE.lo"=quantile(RMSPE, 0.05, na.rm=T), "RMSPE.up"=quantile(RMSPE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T),
              "relRMSPE.est"=mean(relRMSPE, na.rm=T), "relRMSPE.med"=median(relRMSPE, na.rm=T), 
              "relRMSPE.lo"=quantile(relRMSPE, 0.05, na.rm=T), "relRMSPE.up"=quantile(relRMSPE, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSPE.est)  
  
  # -- OPT model
  res_OPT <- res$perf %>%
    data.frame() %>%
    filter(data_ana_model %in% "OPT") %>%
    group_by(data_ana_model, adjust, data_ana_thresh, variable) %>%
    summarise("RMSPE.est"=mean(RMSPE, na.rm=T), "RMSPE.med"=median(RMSPE, na.rm=T), 
              "RMSPE.lo"=quantile(RMSPE, 0.05, na.rm=T), "RMSPE.up"=quantile(RMSPE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T),
              "relRMSPE.est"=mean(relRMSPE, na.rm=T), "relRMSPE.med"=median(relRMSPE, na.rm=T), 
              "relRMSPE.lo"=quantile(relRMSPE, 0.05, na.rm=T), "relRMSPE.up"=quantile(relRMSPE, 0.95, na.rm=T)) %>%
    data.frame() %>%
    arrange(RMSPE.est)   
  
  res_OPT.freq <- res$perf %>%
    filter(data_ana_model %in% "OPT") %>%
    dplyr::count(adjust, data_ana_thresh,  data_ana_t, variable) %>%
    rename(count=n) %>%
    data.frame() %>%
    mutate(perc=round(count/iter*100,2)) %>%
    arrange(variable, data_ana_thresh, desc(count)) 
  
  
  # -- AVG model
  res_AVG <- res$perf %>%
    filter(data_ana_model %in% "AVG") %>%
    select(!data_ana_t) %>%
    group_by(data_ana_model, adjust, data_ana_thresh,  variable) %>%
    summarise("RMSPE.est"=mean(RMSPE, na.rm=T), "RMSPE.med"=median(RMSPE, na.rm=T), 
              "RMSPE.lo"=quantile(RMSPE, 0.05, na.rm=T), "RMSPE.up"=quantile(RMSPE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T),
              "relRMSPE.est"=mean(relRMSPE, na.rm=T), "relRMSPE.med"=median(relRMSPE, na.rm=T), 
              "relRMSPE.lo"=quantile(relRMSPE, 0.05, na.rm=T), "relRMSPE.up"=quantile(relRMSPE, 0.95, na.rm=T)) %>%
    data.frame()
  
  
  # -- FLEX model
  res_FLEX <- res$perf %>%
    filter(data_ana_model %in% "FLEX") %>%
    select(!data_ana_t)  %>%
    group_by(data_ana_model, adjust, data_ana_thresh, variable) %>%
    summarise("RMSPE.est"=mean(RMSPE, na.rm=T), "RMSPE.med"=median(RMSPE, na.rm=T), 
              "RMSPE.lo"=quantile(RMSPE, 0.05, na.rm=T), "RMSPE.up"=quantile(RMSPE, 0.95, na.rm=T),
              "R2.est"=mean(R2, na.rm=T), "R2.med"=median(R2, na.rm=T),
              "R2.lo"=quantile(R2, 0.05, na.rm=T), "R2.up"=quantile(R2, 0.95, na.rm=T),
              "CS.est"=mean(CS, na.rm=T), "CS.med"=median(CS, na.rm=T),
              "CS.lo"=quantile(CS, 0.05, na.rm=T), "CS.up"=quantile(CS, 0.95, na.rm=T),
              "relRMSPE.est"=mean(relRMSPE, na.rm=T), "relRMSPE.med"=median(relRMSPE, na.rm=T), 
              "relRMSPE.lo"=quantile(relRMSPE, 0.05, na.rm=T), "relRMSPE.up"=quantile(relRMSPE, 0.95, na.rm=T)) %>%
    data.frame()
  
  res_FLEX_coeff <- res$more$FLEX_coeff %>%
    group_by(adjust, data_ana_thresh, fda.thresh, variable) %>%
    summarise("coeff.mean.ustd"=mean(fda.est, na.rm=T), 
              "coeff.median.ustd"=median(fda.est, na.rm=T), 
              "coeff.lo.ustd"=quantile(fda.est, 0.05, na.rm=T), 
              "coeff.up.ustd"=quantile(fda.est,0.95, na.rm=T),
              "coeff.mean.std"=mean(fda.est*fda.sd.Xt.pred, na.rm=T), 
              "coeff.median.std"=median(fda.est*fda.sd.Xt.pred, na.rm=T), 
              "coeff.lo.std"=quantile(fda.est*fda.sd.Xt.pred, 0.05, na.rm=T), 
              "coeff.up.std"=quantile(fda.est*fda.sd.Xt.pred,0.95, na.rm=T)) %>%
    data.frame()
  
  res_FLEX_func <- res$more$FLEX_coeff %>%
    group_by(adjust, data_ana_thresh, fda.thresh, variable) %>%
    summarise("aRMSPE.ustd"=calc_rmspe(fda.est, betafun_true),
              "aRMSPE.std"=calc_rmspe(fda.est*fda.sd.Xt.pred, betafun_true)) %>%
    data.frame()
  
  
  # -- Save results 
  tbl_res <- data.frame(rbind(res_oracle, res_null, res_OPT, res_AVG, res_FLEX))
  list_results <- list("scenario"=main.params,
                       "results"=list("tbl_results"= cbind(main.params, tbl_res),
                                      "tbl_OPT_freq"=cbind(main.params, res_OPT.freq),
                                      "tbl_FLEX_coeff"=cbind(main.params, res_FLEX_coeff),
                                      "tbl_FLEX_func"=cbind(main.params, res_FLEX_func)),
                       "iters"=list("res"=cbind(main.params, res$perf),
                                    "fun"=cbind(main.params, res$more$FLEX_coeff)))
  
  
  saveRDS(list_results, here::here(sim.path, paste0(filename , ".rds")))  
}



# ======================== 04. Concatenate simulation results ==================================
evaluate_scenarios <- function(sim.files, sim.path){
  #' Preparation of simulation results for analysis
  
  res <- list()
  for(i in 1:length(sim.files)){
    list.tmp <- readRDS(here::here(sim.path, sim.files[i]))

    res[[i]] <- list("sim"=list.tmp$results$tbl_results, "fun"=list.tmp$results$tbl_FLEX_func, 
                     "tfreq"=list.tmp$results$tbl_OPT_freq, "funcoeff"=list.tmp$results$tbl_FLEX_coeff, 
                     "iters_res"=list.tmp$iters$res, "iters_fun"=list.tmp$iters$fun)
  }
  tbl_scens <- do.call(rbind, lapply(res, function(x) x[[1]]))
  tbl_funs <- do.call(rbind, lapply(res, function(x) x[[2]]))
  tbl_tfreq <- do.call(rbind, lapply(res, function(x) x[[3]]))
  tbl_funcoeff <- do.call(rbind, lapply(res, function(x) x[[4]]))
  tbl_iters_res <- do.call(rbind, lapply(res, function(x) x[[5]]))
  tbl_iters_fun <- do.call(rbind, lapply(res, function(x) x[[6]]))
  
  write.csv(tbl_iters_res,paste0(here::here(sim.path, "/tbl_results_iters.csv")), row.names = FALSE)
  
  out <- list("tbl_scens"=tbl_scens, "tbl_funs"=tbl_funs, "tbl_tfreq"=tbl_tfreq, 
              "tbl_funcoeff"=tbl_funcoeff, 
              "tbl_iters_res"=tbl_iters_res, "tbl_iters_fun"=tbl_iters_fun)
  return(out)
}


# ======================== 05. Report simulation results ==================================
report_simresults <- function(sim.path, filename){
  #' Execute report showcasing simulation results
  
  rmarkdown::render(
    "src/report_main.Rmd",
    output_dir = sim.path,
    output_file = paste0(filename, ".html"))
}
