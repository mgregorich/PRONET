############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Summarize simulation output
############################################



summarize_results <- function(res, params, tseq){
  ltseq <- length(tseq)
  
  res.opt <- data.frame("sim"=1:iter, "opt.thresh"=do.call(rbind, lapply(res, `[[`, 1)))
  
  res.model <- data.frame("sim"=rep(1:params$iter,each=ltseq),do.call(rbind, lapply(res, `[[`, 2)))
  
  res.fitted <- data.frame("sim"=rep(1:params$iter,each=ltseq),do.call(rbind, lapply(res, `[[`, 3)))

  return(paste0("Not yet finished."))
}