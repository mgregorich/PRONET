# ============================================================================ #
# Author: MG
# Date: 07.02.2021
# Info: Simulation study - Setup
# ============================================================================ #


# ======================= Packages & Code ======================================
if(!("pacman" %in% installed.packages())){install.packages("pacman")}
if(!("BiocManager" %in% installed.packages())){install.packages("BiocManager")}
if(!("WGCNA" %in% installed.packages())){BiocManager::install("WGCNA")}
if(!("GO.db" %in% installed.packages())){BiocManager::install("GO.db")}
if(!("looplot" %in% installed.packages())){devtools::install_github("matherealize/looplot")}

pacman::p_load(Rcpp, RcppArmadillo, devtools, MASS, igraph, kableExtra, ggplot2, ggh4x, forcats, 
               gridExtra, ggnewscale, here, stringr, future.apply, parallel, patchwork,
               dplyr,  tidyr, knitr, reshape2, refund, purrr, openxlsx, Hmisc, 
               dtplyr, matrixStats, looplot)
options(dplyr.summarise.inform = FALSE)
sourceCpp(here::here("src","utils.cpp"))


## ======================== Parameters =========================================
#set.seed(666)

# Paths
sim.date <- Sys.Date()
sim.file <- paste0("sim_", sim.date,"/")
sim.path <- here::here("output", sim.file)

source(here::here("src", "functions_main.R"))
source(here::here("src", "functions_aux.R"))


## ========================= Prognostic setting ===================

# -- Data generation
iter = 1                                                                        # iter: number of simulation iterations
p = 116                                                                         # p: number of biomarker nodes
b0 = 10                                                                         # b0: intercept for model
b1 = 20                                                                         # b1: coefficients for network features 
b2 = NA                                                                         # b2: coefficient in multiv. modelling
step_size = 0.01                                                                # step.size: step size for sequential thresholding
data_gen_mech_random_pars = c(0.1, 0.4)
data_gen_mech_universal_par = 0.25

# Varying parameters
n = c(75, 150, 300, 600)                                                        # n: sample size
data_gen_feature = c("cc","cpl")
setting = c("uni")                                                              # setting: univariable or multivariable modelling
data_gen_thresh <- c("weight-based")
data_gen_mech <- c("universal", "random", "flat", "early peak","arc")                                   
epslevel_y = c("none", "moderate", "high")                                        # error term sigma_Y (outcome)
epslevel_g = c("none", "moderate", "high")                                        # error term sigma_G (graph)

scenarios <- expand.grid(
  iter = iter,
  n = n,
  p = p,
  setting = setting,
  data_gen_feature = data_gen_feature,
  data_gen_thresh = data_gen_thresh,
  data_gen_mech = data_gen_mech,
  b0 = b0,
  b1 = b1,
  b2 = b2,
  epslevel_y = epslevel_y,
  epslevel_g = epslevel_g)

# Manual changes to scenarios
scenarios <- scenarios %>%
  arrange(setting, data_gen_feature, data_gen_mech, data_gen_thresh, n) %>%
  mutate(eps_y = 0,
         eps_g = 0) %>%
  mutate(eps_y = case_when(data_gen_mech %in% c("universal", "flat","early peak", "arc") & epslevel_y %in% "moderate" ~ 2.5,
                           data_gen_mech %in% c("universal",  "flat","early peak", "arc") &  epslevel_y %in% "high"~ 5,
                           data_gen_mech %in% c("random") &  epslevel_y %in% "moderate" ~ 3,
                           data_gen_mech %in% c("random") & epslevel_y %in% "high"~ 6,
                           TRUE ~ 0),
         eps_g = case_when(data_gen_mech %in% c("universal", "random",  "flat","early peak", "arc") & epslevel_g %in% "moderate" ~ .15,
                           data_gen_mech %in% c("universal", "random",  "flat","early peak", "arc") & epslevel_g %in% "high" ~ .3,
                           TRUE ~ 0))

print(paste0("Total number of scenarios to be evaluated = ", nrow(scenarios)))

