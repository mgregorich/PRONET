############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data analysis 
############################################


evalLM <- function(x, df){
  # Perform univariable linear regression and extract coeffs, rmse and adjusted r2
  # x="T0.25"; df=data.gvars
  fit.lm <- lm(as.formula(paste0("Y ~ ",x)), data=df)
  sum.lm <- summary(fit.lm)
  rmse <-  sqrt(mean((df$Y-fit.lm$fitted.values)^2))
  ci <- confint(fit.lm)
  
  model.lm <- data.frame("intercept"=as.numeric(fit.lm$coefficients[1]), "int.lo"=ci["(Intercept)",1], "int.upp"= ci["(Intercept)",2],
             "slope"=as.numeric(fit.lm$coefficients[2]), "slope.lo"= ci[x,1], "slope.upp"= ci[x,2],
             "adjr2"=sum.lm$adj.r.squared, "rmse"=rmse)
  out <- list("model"=model.lm, "fitted" = fit.lm$fitted.values)
  return(out)
}


analyse_data <- function(df, tseq){
  # Network data
  data.network <- df[,str_detect(colnames(df), "MI.")]
  
  # CC for threshold sequence
  list.gvars <- lapply(1:nrow(data.network), function(x) evalSSN(eweights=data.network[x,], msize=p, tseq=tseq))
  data.gvars <- data.frame(t(do.call(cbind, list.gvars)))
  colnames(data.gvars) <- paste0("T",seq(0,1,0.05))
  
  # Add outcome Y
  data.gvars <- cbind("Y"=df$Y, data.gvars)
  
  # Pick model with best RMSE
  thresh.names <- colnames(data.gvars)[-1]
  res <- lapply(thresh.names, function(x) evalLM(x, df=data.gvars[,c("Y", x)]))
  res.model <- data.frame("thresh"=thresh.names,do.call(rbind, lapply(res, `[[`, 1)))
  res.fitted <- data.frame("thresh"=thresh.names,do.call(rbind, lapply(res, `[[`, 2)))
  
  opt.thresh <- which.min(res.model$rmse)
  
  out <- list("opt.thresh"=tseq[opt.thresh], "model"=res.model, "y.fit"=res.fitted, "y.true"=df$Y, )
  return(out)
}
# plot(data.sim$FMI, data.gvars$T0.25)  

