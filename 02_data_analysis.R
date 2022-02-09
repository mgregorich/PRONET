# ============================================================================ #
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data analysis 
# ============================================================================ #


analyse_data <- function(df, tseq){
  
  # df=data.test; tseq=thresh.seq
  # Network data
  data.network <- df[,paste0("GE.",1:po)]
  df$Subj <- 1:n
  k = 5
  df$fold <- cvFolds(length(unique(df$Subj)), K=k)$which
  options(dplyr.summarise.inform = FALSE)
  
  
  # CC for threshold sequence
  list.gvars <- lapply(1:nrow(data.network), function(x) data.frame("Subj"=x, wrapperThresholding(eweights=data.network[x,], msize=p, tseq=tseq)))
  data.gvars <- data.frame((do.call(rbind, list.gvars)))

  # Add outcome Y
  data.gvars <- merge(df[,c("Subj","Y", "fold")], data.gvars, by="Subj") %>%
    mutate(Value=ifelse(is.nan(Value), NA, Value)) %>%
    filter(!Variable %in% "ncon")

  
  # ------ Pick model with best RMSE
  data.bRMSE <- data.gvars  %>%
    group_by(SparsMethod, ThreshMethod, Thresh, Variable) %>%
    mutate("Y.pred"=evalLM(Y, X=Value, fold=fold, k=k)) %>%
    mutate(RMSE=sqrt(mean((Y-Y.pred)^2))) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    slice(which.min(RMSE)) %>%
    select(!c(Y,Y.pred, Subj)) %>%
    mutate("AnaMethod"="bRMSE")
  
  # ------ Average feature across threshold sequence
  data.AVG <- data.gvars %>%
    group_by(SparsMethod, ThreshMethod, Variable, Subj, Y, fold) %>%
    summarise("Value.avg"=mean(Value, na.rm=T)) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    mutate("Y.pred"=evalLM(Y, X=Value.avg, fold=fold, k=k)) %>%
    summarise(RMSE=sqrt(mean((Y-Y.pred)^2))) %>%
    mutate("AnaMethod"="AVG",
           "Thresh"=NA)
  
  # ------- Functional data analysis approach
  data.FDA <- data.gvars %>%
    pivot_wider(values_from = Value, names_from = Thresh) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    do("RMSE"=evalPFR(x= ., fold=fold, k=k, tseq=tseq)) %>%
    mutate(RMSE=unlist(RMSE),Thresh=NA,
           "AnaMethod"="FDA")
  
  # ------- Results
  out <- data.frame(rbind(data.bRMSE[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE")],
                   data.AVG[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE")],
                   data.FDA[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE")]))
  return(out)
}

