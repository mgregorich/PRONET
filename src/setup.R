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

pacman::p_load(mvtnorm, igraph, NetworkToolbox, Rcpp, RcppArmadillo, RcppEigen, MASS, lqmm, 
               kableExtra, ggplot2, ggh4x, forcats, gridExtra, here,
               stringr, future.apply, parallel, dplyr, tidyr, knitr, reshape2,
               refund, broom, cvTools, concreg, fda, purrr, openxlsx, DT,
               dtplyr, profvis, matrixStats, Rfast, looplot)
options(dplyr.summarise.inform = FALSE)
sourceCpp(here::here("src","thresholding.cpp"))

## ======================== Parameters =========================================
set.seed(666)

# -- Data generation
iter = 10                                                                      # number of simulation iterations
q = 2                                                                           # q: number of covariates; 
b0 = 10                                                                         # intercept for model
b1 = 10                                                                         # coefficients for network features 

# Varying parameters
n = c(125, 250, 500)                                                            # n: sample size
p = c(50, 100, 200)                                                             # p: number of biomarker nodes
dg.thresh = list("single"=c(0.25),                                              # Sparsification threshold for data gen
                 "random"=seq(0.1,0.4,0.02),
                 "func"="flat",
                 "func"="half-sine",
                 "func"="sine")                                   
eps.y = c(0, .5,  1)                                                            # error term sigma_Y (outcome)
eps.g = c(0, .05, .1)                                                           # error term sigma_G (graph)

# -- Parameter distribution for edge weights ~ beta(a,b)
beta.params = list(c(2,5))                                                      # shape params of beta distribution
alpha0.params = list("norm"=c("mean"=5, "sd"=2.5))                              # stat params of normal distributed alpha0
alpha12.params = list("unif"=c("min"=0, "max"=2))                               # stat params of uniform distributed alpha1 and alpha2
X1.params = list("norm"=c("mean"=0, "sd"=2))                                    # stat params of normal distributed latent processes X1 and X2
X2.params = list("binom"=0.5)                                                   # stat params of normal distributed latent processes X1 and X2
excel = F                                                                       # generate additional excel file with scen results

scenarios <- expand.grid(
  iter = iter,
  n = n,
  q = q,
  p = p,
  dg.thresh = dg.thresh,
  beta.params = beta.params,
  alpha0.params = alpha0.params,
  alpha12.params = alpha12.params,
  X1.params = X1.params,
  X2.params = X2.params,
  b0 = b0,
  b1 = b1,
  eps.y = eps.y,
  eps.g = eps.g,
  excel = excel 
)

print(paste0("Total number of scenarios to be evaluated = ", nrow(scenarios)))


