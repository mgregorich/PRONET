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

pacman::p_load(mvtnorm, Rcpp, MASS, Matrix, igraph,
               kableExtra, ggplot2, ggh4x, forcats, gridExtra, ggnewscale, here,
               stringr, future.apply, parallel, dplyr,  tidyr, knitr, reshape2,
               refund, broom, cvTools, concreg,  purrr, openxlsx,
               dtplyr, profvis, matrixStats, looplot, corpcor)
options(dplyr.summarise.inform = FALSE)
sourceCpp(here::here("src","utils.cpp"))



## ======================== Parameters =========================================
set.seed(666)

# Paths
sim.date <- Sys.Date()
sim.file <- paste0("sim_", sim.date,"/")
sim.path <- here::here("output", sim.file)

if(!dir.exists(here::here("output"))){dir.create(here::here("output"))}
if(dir.exists(sim.path)){invisible(do.call(file.remove, list(list.files(sim.path, full.names = T))))
}else{dir.create(sim.path)}

source(here::here("src", "functions_main.R"))
source(here::here("src", "functions_aux.R"))

# -- Data generation
iter = 3                                                                      # number of simulation iterations
p = 150                                                                         # p: number of biomarker nodes
q = 2                                                                           # q: number of covariates; 
b0 = 10                                                                         # intercept for model
b1 = 10                                                                   # coefficients for network features 
b2 = NA
step.size = 0.01 

# Varying parameters
n = c(75, 150, 300)                                                             # n: sample size
setting = c("uni", "latent", "multi")
network.model = c("random")
dg.spars = c("weight-based")
dg.thresh = list("single"=c(0.25),                                              # Sparsification threshold for data gen
                 "random"=c(0.1,0.4),
                 "flat"="flat",
                 "half-sine"="half-sine",
                 "sine"="sine")                                   
epslevel.y = c("none", "medium", "high")                                        # error term sigma_Y (outcome)
epslevel.g = c("none", "medium", "high")                                        # error term sigma_G (graph)

# -- Parameter distribution for edge weights ~ beta(a,b)
beta.params = list(c(2,6))                                                      # shape params of beta distribution
alpha0.params = list("norm"=c("mean"=5, "sd"=2.5))                              # stat params of normal distributed alpha0
alpha12.params = list("unif"=c("min"=0, "max"=2))                               # stat params of uniform distributed alpha1 and alpha2
Z1.params = list("norm"=c("mean"=0, "sd"=2))                                    # stat params of normal distributed latent processes Z1 and Z2
Z2.params = list("binom"=0.5)                                                   # stat params of normal distributed latent processes Z1 and Z2
excel = F                                                                       # generate additional excel file with scen results

scenarios <- expand.grid(
  setting = setting,
  iter = iter,
  n = n,
  q = q,
  p = p,
  network.model=network.model,
  dg.spars = dg.spars,
  dg.thresh = dg.thresh,
  beta.params = beta.params,
  alpha0.params = alpha0.params,
  alpha12.params = alpha12.params,
  Z1.params = Z1.params,
  Z2.params = Z2.params,
  b0 = b0,
  b1 = b1,
  b2 = b2,
  epslevel.y = epslevel.y,
  epslevel.g = epslevel.g,
  step.size = step.size,
  excel = excel 
)

print(paste0("Total number of scenarios to be evaluated = ", nrow(scenarios)))

# 
scenarios <- scenarios %>%
  arrange(setting,n) %>%
  mutate(eps.y = 0,
         eps.g = 0) %>%
  mutate(eps.g = case_when(names(dg.thresh) %in% c("single", "flat", "sine", "half-sine") & epslevel.g %in% "medium"~ 1,
                           names(dg.thresh) %in% c("single", "flat", "sine", "half-sine") & epslevel.g %in% "high" ~ 1.5,
                           names(dg.thresh) %in% c("random") & epslevel.g %in% "medium" ~ 1.5,
                           names(dg.thresh) %in% c("random") & epslevel.g %in% "high" ~ 2.5,
                           TRUE ~ 0),
         eps.y = case_when(names(dg.thresh) %in% c("random", "sine") & epslevel.y %in% "medium" ~ 2,
                           names(dg.thresh) %in% c("random", "sine") & epslevel.y %in% "high"~ 4,
                           names(dg.thresh) %in% c("flat") & epslevel.y %in% "medium" ~ 1,
                           names(dg.thresh) %in% c("flat") & epslevel.y %in% "high"~ 2.5,
                           names(dg.thresh) %in% c("single", "half-sine") & epslevel.y %in% "medium"~ 1.5,
                           names(dg.thresh) %in% c("single", "half-sine") & epslevel.y %in% "high"~ 3,
                           TRUE ~ 0),
         eps.g = case_when(names(dg.thresh) %in% c("half-sine","single", "sine", "flat") & setting %in% "latent" ~ eps.g*1.5,
                           names(dg.thresh) %in% c("random") & setting %in% "latent" ~ eps.g*1,
                           TRUE ~ eps.g),
         eps.y = case_when(names(dg.thresh) %in% c("random", "half-sine","single", "sine", "flat") & setting %in% "latent" ~ eps.y*1.25,
                           TRUE ~ eps.y),
         eps.g = case_when(names(dg.thresh) %in% c("random", "half-sine","single", "sine", "flat") & setting %in% "multi" ~ eps.g*1,
                           TRUE ~ eps.g),
         eps.y = case_when(names(dg.thresh) %in% c("random", "single", "sine") & setting %in% "multi" ~ eps.y*1,
                           names(dg.thresh) %in% c("flat", "half-sine") & setting %in% "multi" ~ eps.y*1,
                           TRUE ~ eps.y))
# scenarios[scenarios$dg.spars=="density-based" & names(scenarios$dg.thresh) %in% "single",]$dg.thresh <- 0.75
# scenarios[scenarios$dg.spars=="density-based" & names(scenarios$dg.thresh) %in% "random",]$dg.thresh <- lapply(scenarios[scenarios$dg.spars=="density-based" & names(scenarios$dg.thresh) %in% "random",]$dg.thresh, function(x) x<-c(0.6,0.9))
#                            



