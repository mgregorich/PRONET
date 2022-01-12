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
    mutate(Value=ifelse(is.nan(Value), NA, Value))
  
  
  # Pick model with best RMSE
  data.gvars <- data.gvars %>%
    group_by(SparsMethod, ThreshMethod, Variable, Thresh) %>%
    mutate("Y.pred"=evalLM(Y, Value)) %>%
    mutate(rmse=sqrt(mean((Y-Y.pred)^2)))
  
  opt.thresh <- which.min(data.gvars$rmse)
  
  out <- data.gvars[opt.thresh,c("SparsMethod", "ThreshMethod", "Thresh", "Variable","rmse")]
  return(out)
}
# plot(data.sim$FMI, data.gvars$T0.25)  

