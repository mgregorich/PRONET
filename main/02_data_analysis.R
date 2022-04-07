# ============================================================================ #
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data analysis 
# ============================================================================ #


analyse_data <- function(df, tseq, true.params, k=5){
  
  # df=data.test; tseq=da.thresh; k=5; true.params=list("SparsMethod"="weight-based", "ThreshMethod"="trim", "Thresh"=sparams$sthresh)
  # Network data
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

  # ------ Oracle model
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
  
  
  # ------ Pick model with best RMSE
  data.bRMSE <- data.gvars  %>% 
    tibble() %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM_dCV(data.lm=df,k=k))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="bRMSE")
  
  
  # ------ Average feature across threshold sequence
  data.AVG <- data.gvars %>%
    group_by(SparsMethod, ThreshMethod, Variable, Subj, Y, fold) %>%
    summarise("Value.avg"=mean(Value, na.rm=T)) %>%
    rename(Value=Value.avg) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    nest() %>%
    mutate(res=lapply(data, function(df) evalLM(data.lm=df, k=k))) %>%
    unnest(res, keep_empty = T) %>%
    mutate("AnaMethod"="AVG")
  
  
  # ------- Functional data analysis approach
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
  
  # ------- Results
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

