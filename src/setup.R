# ============================================================================ #
# Author: MG
# Date: 07.02.2021
# Info: Simulation study - Setup
# ============================================================================ #


# ======================= Packages & Code ======================================
if(!("BiocManager" %in% installed.packages())){install.packages("BiocManager")}
if(!("WGCNA" %in% installed.packages())){BiocManager::install("WGCNA")}
if(!("GO.db" %in% installed.packages())){BiocManager::install("GO.db")}

pacman::p_load(mvtnorm, igraph, NetworkToolbox, Rcpp, RcppEigen, MASS, lqmm, 
               kableExtra, ggplot2, ggh4x, forcats, gridExtra, here,
               stringr, future.apply, parallel, dplyr, tidyr, knitr, reshape2,
               refund, broom, cvTools, concreg, fda, purrr, openxlsx, DT,
               dtplyr, profvis, matrixStats, Rfast)
options(dplyr.summarise.inform = FALSE)

## ======================== Parameters =========================================
set.seed(666)

# -- Data generation
iter = 5                                                                        # number of simulation iterations
n = c(100, 250, 500)                                                                 # n: sample size
q = 2                                                                           # q: number of covariates; 
p = c(25, 50, 100)                                                                  # p: number of biomarker nodes
dg.thresh = list("single"=c(0.25),                                             # Sparsification threshold for data gen
                 "random"=seq(0.1,0.4,0.05),
                 "func"="quadratic")                                   
b0 = 10                                                                         # intercept for model
b1 = 10                                                                         # coefficients for network features 
eps.y = c(0, .5,  1, 1.5)                                                       # error term sigma_Y (outcome)
eps.g = c(0, .05, .1)                                                           # error term sigma_G (graph)
excel = F                                                                       # generate additional excel file with scen results

# -- Parameter distribution for edge weights ~ beta(a,b)
beta.params = list(c(2,5))                                                      # c(4,2),c(4,4)  # shape params of beta distribution
alpha0.params = list("norm"=c("mean"=5, "sd"=2.5))                              # stat params of normal distributed alpha0
alpha12.params = list("unif"=c("min"=0, "max"=2))                               # stat params of uniform distributed alpha1 and alpha2
X1.params = list("norm"=c("mean"=0, "sd"=2))                                    # stat params of normal distributed latent processes X1 and X2
X2.params = list("binom"=0.5)                                                   # stat params of normal distributed latent processes X1 and X2

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

# Barabasi-Albert model with quadratic preferential attachment for Bernoulli graph
BA.graph <- sample_pa(n=p, power=1, m=30, directed = F)                         # increase m to increase density
plot(BA.graph, vertex.label.dist=1.5, vertex.size=3, vertex.label=NA)
edge_density(BA.graph)

