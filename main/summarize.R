#===============================================================================#
# Author: MG
# Date: 21.04.2021
# Info: Summarize and evaluate all executed scenarios
#===============================================================================#



# ======================== Set up ==============================================
rm(list=ls())

# Functions
source("main/functions_main.R")
source("main/functions_aux.R")
source("main/setup.R")


# Paths
out.path <- "output/"
sim.file <- "sim_2022-04-21/"              #paste0("sim_", Sys.Date(),"/")
sim.path <- paste0(out.path, sim.file)


# Summarize all scenarios
sim.files <- list.files(sim.path, pattern = "sim_")
sim.all <- lapply(sim.files, function(x) cbind_results(x, sim.path))
tbl_scens <- do.call(rbind, sim.all)

saveRDS(tbl_scens, paste0(sim.path, "tbl_scenarios_results.rds"))  
write.xlsx(tbl_scens, paste0(sim.path, "tbl_scenarios_results.xlsx"), overwrite = T)

# Preprocess
tbl_scens <- tbl_scens %>%
  mutate_at(vars(dg.thresh), as.factor)

tbl_scens %>%
  filter(ThreshMethod %in% "trim" & eps.y == 0.5 & eps.g==0.025 &!(AnaMethod %in% "Oracle")) %>%
  ggplot(aes(x=dg.thresh, y=RMSE.est, group=AnaMethod, col=AnaMethod)) +
  geom_ribbon(aes(x=dg.thresh, ymax=RMSE.up, ymin=RMSE.lo, fill=AnaMethod), linetype=2,alpha=.05) +
  geom_line() +
  geom_point() +
  scale_x_discrete("True threshold in data generation") +
  scale_y_continuous("RMSE") +
  facet_grid(SparsMethod~beta.params) +
  theme_bw()

tbl_scens %>%
  filter(ThreshMethod %in% "trim" & eps.y == 0.75  & !(AnaMethod %in% "Oracle")) %>%
  ggplot(aes(x=dg.thresh, y=RMSE.est, group=AnaMethod, col=AnaMethod)) +
  geom_ribbon(aes(x=dg.thresh, ymax=RMSE.up, ymin=RMSE.lo, fill=AnaMethod), linetype=2,alpha=.05) +
  geom_line() +
  scale_x_discrete("True threshold in data generation") +
  scale_y_continuous("RMSE") +
  geom_point() +
  facet_grid(SparsMethod~beta.params) +
  theme_bw()

tbl_scens %>%
  filter(ThreshMethod %in% "trim" & eps.y ==1 & eps.g==0.025 & !(AnaMethod %in% "Oracle")) %>%
  ggplot(aes(x=dg.thresh, y=RMSE.est, group=AnaMethod, col=AnaMethod)) +
  geom_ribbon(aes(x=dg.thresh, ymax=RMSE.up, ymin=RMSE.lo, fill=AnaMethod), linetype=2,alpha=.05) +
  geom_line() +
  geom_point() +
  scale_x_discrete("True threshold in data generation") +
  scale_y_continuous("RMSE") +
  facet_grid(SparsMethod~beta.params) +
  theme_bw()

#-----------------

tbl_scens %>%
  filter(ThreshMethod %in% "trim" & dg.thresh ==0.4 & eps.g==0.025 & !(AnaMethod %in% "Oracle")) %>%
  ggplot(aes(x=eps.y, y=RMSE.est, group=AnaMethod, col=AnaMethod)) +
  geom_ribbon(aes(x=eps.y, ymax=RMSE.up, ymin=RMSE.lo, fill=AnaMethod), linetype=2,alpha=.05) +
  geom_line() +
  geom_point() +
  scale_x_discrete(expression(epsilon[Y])) +
  scale_y_continuous("RMSE") +
  facet_grid(SparsMethod~beta.params) +
  theme_bw()
