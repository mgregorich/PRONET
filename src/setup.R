# ============================================================================ #
# Author: MG
# Date: 07.02.2021
# Info: Simulation study - Setup
# ============================================================================ #


# ======================= Packages & Code ======================================
pacman::p_load(mvtnorm, igraph, NetworkToolbox, Rcpp, RcppEigen, MASS, lqmm, 
               kableExtra, ggplot2, ggh4x, forcats, gridExtra, here,
               stringr, future.apply, parallel, dplyr, tidyr, knitr, reshape2,
               refund, refund.shiny, broom, cvTools, concreg, fda, purrr, openxlsx, DT,
               dtplyr, profvis)
options(dplyr.summarise.inform = FALSE)

## ======================== Parameters =========================================
set.seed(666)
    
# -- Data generation
iter = 2
5                                                                       # number of simulation iterations
n = 250                                                                         # n: sample size
q = 2                                                                           # q: number of covariates; 
p = 50                                                                          # p: number of biomarker nodes
dg.thresh = c(0,.25,.5,.75)                                                      # Sparsification threshold for data gen
da.thresh = list(seq(0.1,0.5,.02))                                                # Sparsification sequence for data ana
b0 = 10                                                                         # intercept for model
b1 = 5                                                                        # coefficients for network features 
eps.y = c(0, .5, 1, 1.5)                                                        # error term sigma_Y (outcome)
eps.g = c(0, .05, .1, .15, .2)                                                      # error term sigma_G (graph)
report = F                                                                      # generate report for scenario
excel = F                                                                       # generate additional excel file with scen results

# -- Parameter distribution for edge weights ~ beta(a,b)
beta.params = list(c(4,2),c(4,4),c(2,4))                                        # shape params of beta distribution
alpha0.params = list("norm"=c("mean"=3, "sd"=1))                                # stat params of normal distributed alpha0
alpha12.params = list("unif"=c("min"=0, "max"=2))                               # stat params of uniform distributed alpha1 and alpha2
X.params = list("norm"=c("mean"=0, "sd"=1))                                     # stat params of normal distributed latent processes X1 and X2

scenarios <- expand.grid(
  iter = iter,
  n = n,
  q = q,
  p = p,
  dg.thresh = dg.thresh,
  da.thresh = da.thresh,
  beta.params = beta.params,
  alpha0.params = alpha0.params,
  alpha12.params = alpha12.params,
  X.params = X.params,
  b0 = b0,
  b1 = b1,
  eps.y = eps.y,
  eps.g = eps.g,
  report = report,
  excel = excel 
)

print(paste0("Total number of scenarios to be evaluated = ", nrow(scenarios)))

