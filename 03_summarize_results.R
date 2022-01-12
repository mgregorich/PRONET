############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Summarize simulation output
############################################



summarize_results <- function(res, sparams, tseq){
  # res=results.sim; sparams=sparams; tseq=thresh.seq
  
    res <- res %>%
    count(SparsMethod, ThreshMethod, Thresh, Variable) %>%
    rename(count=n) %>%
    data.frame() %>%
    mutate(perc=count/sum(count)*100) %>%
    arrange(desc(count))
  
  out <- list("tbl_summary"=res)
  return(out)
}
