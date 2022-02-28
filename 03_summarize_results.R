# ============================================================================ #
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Summarize simulation output
# ============================================================================ #


res <- results.sim %>%
  mutate(iter=rep(1:iter, each=18))

# ------------- Best RMSE ----------------

res.bRMSE <- res %>%
    filter(AnaMethod %in% "bRMSE") %>%
    group_by(SparsMethod, ThreshMethod, Thresh, Variable) %>%
    summarise("RMSE.mean"=mean(RMSE), "R2.mean"=mean(R2), "CS.mean"=mean(CS)) %>%
    data.frame() %>%
    arrange(RMSE.mean) %>%
    filter(row_number() %in% 1:3)  

res.bRMSE.freq <- res %>%
  filter(AnaMethod %in% "bRMSE") %>%
  group_by(Variable, iter) %>%
  slice_min(RMSE) %>%
  ungroup() %>%
  select(!iter) %>%
  count(SparsMethod, ThreshMethod, Thresh, Variable) %>%
  rename(count=n) %>%
  data.frame() %>%
  mutate(perc=count/iter*100) %>%
  group_by(SparsMethod, Thresh, Variable, count) %>%
  slice_head() %>%
  arrange(Variable, desc(count)) 

  
# --------- Averaging over threshold sequence

res.AVG <- res %>%
  filter(AnaMethod %in% "AVG") %>%
  select(!Thresh) %>%
  group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
  summarise("RMSE.mean"=mean(RMSE), "R2.mean"=mean(R2), "CS.mean"=mean(CS)) %>%
  pivot_wider(values_from = RMSE.mean, names_from = ThreshMethod) %>%
  arrange(Variable) %>%
  relocate(bin, .after=Variable)
  

# ----------- FDA ----------------
res.FDA <- res %>%
  filter(AnaMethod %in% "FDA") %>%
  select(!Thresh)  %>%
  group_by(AnaMethod, SparsMethod, ThreshMethod, Variable) %>%
  summarise("RMSE.mean"=mean(RMSE), "R2.mean"=mean(R2), "CS.mean"=mean(CS)) %>%
  pivot_wider(values_from = RMSE.mean, names_from = ThreshMethod) %>%
  arrange(Variable) %>%
  relocate(bin, .after=Variable)

list_results <- list("tbl_bRMSE"=res.bRMSE, "tbl_bRMSE_freq"=res.bRMSE.freq, "tbl_AVG"=res.AVG, "tbl_FDA"=res.FDA)
saveRDS(list_results, paste0(sim.path, "Robject_results.rds"))
