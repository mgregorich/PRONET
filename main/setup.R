# ============================================================================ #
# Author: MG
# Date: 07.02.2021
# Info: Simulation study - Setup
# ============================================================================ #


# ======================= Packages & Code ======================================
pacman::p_load(mvtnorm, igraph, NetworkToolbox, Rcpp, RcppEigen, MASS, lqmm, 
               kableExtra, ggplot2, forcats, gridExtra, here,
               stringr, future.apply, parallel, dplyr, tidyr, knitr, reshape2,
               refund, refund.shiny, broom, cvTools, concreg, fda, purrr)

## ======================== Parameters =========================================
set.seed(666)

# -- Data generation
iter = 10                                                                       # number of simulation iterations
n = 250                                                                         # n: sample size
q = 2                                                                           # q: number of covariates; 
p = 50                                                                          # p: number of biomarker nodes
dg.thresh = 0.35                                                                # Sparsification threshold for data gen
da.thresh = seq(0,1,0.05)                                                       # Sparsification sequence for data ana
po = (p-1)*p/2                                                                  # po: number of possible undirected edges
beta0 = 10                                                                      # intercept for model
xbeta = NA                                                                      # coefficients for covariate X
gbeta = 5                                                                       # coefficients for network features fmi
eps.y = .25                                                                     # error term sigma_Y (outcome)
eps.g = .025                                                                    # error term sigma_G (graph)

beta.par1 = 4
beta.par2 = 2

# -- Barabasi-Albert model for Bernoulli graph
BA.graph <- sample_pa(n=p, power=2, m=20, directed = F)                         # increase m to increase density

# -- Parameter distribution for edge weights ~ beta(a,b)
distr.params=list("beta"=c("shape1"=beta.par1, "shape2"=beta.par2), 
                  "alpha0.norm"=c("mean"=3, "sd"=1),
                  "alpha12.unif"=c("min"=0, "max"=2), 
                  "X.norm"=c("mean"=0, "sd"=1)) 
true.params = list("SparsMethod"="weight-based",
                   "ThreshMethod"="trim",
                   "Thresh"=dg.thresh)

# -- Save all relevant parameters in list
main.params <- list(iter=iter,n=n, p=p, q=q,
                beta0=beta0, xbeta=xbeta, gbeta=gbeta, sthresh=dg.thresh, 
                eps.y=eps.y, eps.g=eps.g)
