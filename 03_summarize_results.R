############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Summarize simulation output
############################################



summarize_results <- function(res, sparams, tseq){
  # res=results.sim; sparams=sparams; tseq=thresh.seq
  ltseq <- length(tseq)
  
  # Opt thresh
  res.opt <- data.frame("sim"=1:iter, "opt.thresh"=do.call(rbind, lapply(res, `[[`, 1)))
  tmp1 <- res.opt %>%
    summarize( mean.thresh=mean(opt.thresh),
               n.corr = sum(res.opt$opt.thresh == sparams$thresh),
               perc.corr = sum(res.opt$opt.thresh == sparams$thresh)/iter*100) %>%
    mutate_all(round,2)
  
  
  # Coef res
  res.model <- data.frame("sim"=rep(1:sparams$iter,each=ltseq),do.call(rbind, lapply(res, `[[`, 2)))
  tmp2 <- res.model %>% 
    group_by(thresh) %>%
    summarize(icpt_mean.hat = mean(intercept, na.rm = T),
              icpt_bias = icpt_mean.hat - sparams$beta0,
              icpt_relbias = ((icpt_mean.hat - sparams$beta0)/sparams$beta0)*100,
              icpt_var = var(intercept,na.rm=T),
              icpt_rmse = sqrt(mean((intercept - sparams$beta0) ^ 2, na.rm=T)),
              coef_mean.hat = mean(slope, na.rm = T),
              coef_bias = coef_mean.hat - sparams$gbeta,
              coef_relbias = ((coef_mean.hat - sparams$gbeta)/sparams$gbeta)*100,
              coef_var = var(slope,na.rm=T),
              coef_rmse = sqrt(mean((slope - sparams$gbeta) ^ 2, na.rm=T))) %>% 
    mutate_at(2:11, round, 3) 
  
  
  res.fitted <- data.frame("sim"=rep(1:sparams$iter,each=ltseq),do.call(rbind, lapply(res, `[[`, 3)))
  res.true <- data.frame("sim"=rep(1:sparams$iter,each=ltseq),do.call(rbind, lapply(res, `[[`, 4)))
                
  out <- list("tbl_thresh"=tmp1, "tbl_coef"=tmp2)
  return(out)
}
