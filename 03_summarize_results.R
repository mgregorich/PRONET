# ============================================================================ #
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Summarize simulation output
# ============================================================================ #



# ========================== Extraction of results =============================
res <- list()
# Overall performance comparison between methods 
res$perf <- do.call(rbind, lapply(results.sim,  function(x) x[[1]])) %>%
  mutate(iter=rep(1:iter, each=19))

res$more <- list()
# Coefficient function of the functional data approach
res$more$FDA.coeff <- do.call(rbind, lapply(1:length(results.sim),  function(x) data.frame(iter=x, results.sim[[x]]$more$FDA.coeff)))


# ================================= Methods ====================================

# ----------- Oracle ------------------
res.oracle <- res$perf %>%
  data.frame() %>%
  filter(AnaMethod %in% "Oracle") %>%
  group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
  summarise("RMSE.est"=mean(RMSE), "RMSE.lo"=quantile(RMSE, 0.05), "RMSE.up"=quantile(RMSE, 0.95),
            "R2.est"=mean(R2),"R2.lo"=quantile(R2, 0.05), "R2.up"=quantile(R2, 0.95),
            "CS.est"=mean(CS),"CS.lo"=quantile(CS, 0.05), "CS.up"=quantile(CS, 0.95)) %>%
  data.frame() %>%
  arrange(RMSE.est)  


# ----------- Best RMSE --------------
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

  
# --------- Averaging over threshold sequence ----------
res.AVG <- res$perf %>%
  filter(AnaMethod %in% "AVG") %>%
  select(!Thresh) %>%
  group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
  summarise("RMSE.est"=mean(RMSE), "RMSE.lo"=quantile(RMSE, 0.05), "RMSE.up"=quantile(RMSE, 0.95),
            "R2.est"=mean(R2),"R2.lo"=quantile(R2, 0.05), "R2.up"=quantile(R2, 0.95),
            "CS.est"=mean(CS),"CS.lo"=quantile(CS, 0.05), "CS.up"=quantile(CS, 0.95)) 
  

# ----------- FDA ----------------
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


# ========================== Save results ======================================
list_results <- list("tbl_oracle"=res.oracle,
                     "tbl_bRMSE"=res.bRMSE, "tbl_bRMSE_freq"=res.bRMSE.freq, 
                     "tbl_AVG"=res.AVG, 
                     "tbl_FDA"=res.FDA, "tbl_FDA_coeff"=res.FDA.coeff,
                     "data_example" =results.sim[[1]]$data)
saveRDS(list_results, paste0(sim.path, "Robject_results.rds"))

