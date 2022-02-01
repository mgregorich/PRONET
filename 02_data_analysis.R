############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data analysis 
############################################


evalLM <- function(Y, X){
  # Perform univariable linear regression and extract coeffs, rmse and adjusted r2
  # df<- data.gvars %>% filter(SparsMethod == "density-based"& ThreshMethod == "trim"& Variable == "cpl"& Thresh == 0); Y=df$Y; X=df$Value
  
  if(sum(X, na.rm=T)==0){
    out <- rep(NA, length(X))
  }else{
    fit.lm <- lm(Y~X, na.action = "na.exclude")
    # sum.lm <- summary(fit.lm)
    # rmse <-  sqrt(mean((df$Y-fit.lm$fitted.values)^2))
    # ci <- confint(fit.lm)
    # 
    # model.lm <- data.frame("intercept"=as.numeric(fit.lm$coefficients[1]), "int.lo"=ci["(Intercept)",1], "int.upp"= ci["(Intercept)",2],
    #            "slope"=as.numeric(fit.lm$coefficients[2]), "slope.lo"= ci[x,1], "slope.upp"= ci[x,2],
    #            "adjr2"=sum.lm$adj.r.squared, "rmse"=rmse)
    out <- fitted(fit.lm)
  }

  return(out)
}

evalPFR <- function(Y,X, tseq){
 # X=data.wider[data.wider$SparsMethod %in% "density-based" & data.wider$ThreshMethod %in% "trim" & data.wider$Variable %in% "cc",]; Y=Y.obs
  
  flength <- length(tseq)
  X <- as.matrix.data.frame(X[,(ncol(X)-flength+1):ncol(X)])
  fit.fda <- refund::pfr(Y ~ lf(X, k = 5, presmooth="bspline"), 
                         family="gaussian")
  fit.fda.sum <- summary(fit.fda)
  y_pred = c(predict(fit.fda, newdata=list(X=X), type="response"))
  RMSE=sqrt(mean((Y-y_pred)^2))
  return(RMSE)
}


analyse_data <- function(df, tseq){
  # df=data.test; tseq=thresh.seq
  # Network data
  data.network <- df[,paste0("GE.",1:po)]
  df$Subj <- 1:n
  
  # CC for threshold sequence
  list.gvars <- lapply(1:nrow(data.network), function(x) data.frame("Subj"=x, wrapperThresholding(eweights=data.network[x,], msize=p, tseq=tseq)))
  data.gvars <- data.frame((do.call(rbind, list.gvars)))

  # Add outcome Y
  data.gvars <- merge(df[,c("Subj","Y")], data.gvars, by="Subj") %>%
    mutate(Value=ifelse(is.nan(Value), NA, Value)) %>%
    filter(!Variable %in% "ncon")
  
  # ------ Pick model with best RMSE
  data.bRMSE <- data.gvars  %>%
    group_by(SparsMethod, ThreshMethod, Variable, Thresh) %>%
    mutate("Y.pred"=evalLM(Y, Value)) %>%
    mutate(RMSE=sqrt(mean((Y-Y.pred)^2))) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    slice(which.min(RMSE)) %>%
    select(!c(Y,Y.pred, Subj)) %>%
    mutate("AnaMethod"="bRMSE")
  
  # ------ Average feature across threshold sequence
  data.AVG <- data.gvars %>%
    group_by(SparsMethod, ThreshMethod, Variable, Subj,Y) %>%
    summarise("Value.avg"=mean(Value, na.rm=T)) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    mutate("Y.pred"=evalLM(Y, Value.avg)) %>%
    summarise(RMSE=sqrt(mean((Y-Y.pred)^2))) %>%
    mutate("AnaMethod"="AVG",
           "Thresh"=NA)
  
  # ------- Functional data analysis approach
  Y.obs <- data.gvars[!duplicated(data.gvars[, c("Subj", "Y")]),]$Y
  data.FDA <- data.gvars %>%
    pivot_wider(values_from = Value, names_from = Thresh) %>%
    group_by(SparsMethod, ThreshMethod, Variable) %>%
    do("RMSE"=evalPFR(Y.obs, X= ., tseq=tseq)) %>%
    mutate(RMSE=unlist(RMSE),Thresh=NA,
           "AnaMethod"="FDA")
  
  # ------- Results
  out <- data.frame(rbind(data.bRMSE[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE")],
                   data.AVG[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE")],
                   data.FDA[,c("AnaMethod","SparsMethod", "ThreshMethod", "Thresh","Variable","RMSE")]))
  return(out)
}

