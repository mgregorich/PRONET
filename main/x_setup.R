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

# -- Parameter distribution for edge weights ~ beta(a,b)
distr.params=list("beta"=c("shape1"=beta.par1, "shape2"=beta.par2), 
                  "alpha0.norm"=c("mean"=3, "sd"=1),
                  "alpha12.unif"=c("min"=0, "max"=2), 
                  "X.norm"=c("mean"=0, "sd"=1)) 
true.params = list("SparsMethod"="weight-based",
                   "ThreshMethod"="trim",
                   "Thresh"=dg.thresh)


# ======================  Network setup ========================================
# -- Barabasi-Albert model for Bernoulli graph
BA.graph <- sample_pa(n=p, power=2, m=20, directed = F)                         # increase m to increase density
BA.strc <- as.matrix(as_adjacency_matrix(BA.graph))
# edge_density(BA.graph)

# -- Edge weights ~ beta(a,b)
eta.params <- calc_eta_mean_and_var(alpha0.norm.pars=distr.params$alpha0.norm, 
                                    X.norm.pars=distr.params$X.norm,
                                    alpha12.unif.pars=distr.params$alpha12.unif)
alpha0 <- BA.strc[lower.tri(BA.strc)]
alpha0[alpha0==1] <- rnorm(sum(alpha0), distr.params$alpha0.norm[1], 
                           distr.params$alpha0.norm[2])
omega.imat=matrix(alpha0,1,po, byrow = T)

alpha12 <- rep(BA.strc[lower.tri(BA.strc)],q)
alpha12[alpha12==1] <- runif(sum(alpha12), distr.params$alpha12.unif[1], 
                             distr.params$alpha12.unif[2])                             
alpha12.wmat=matrix(alpha12,q,po, byrow = T)
alpha=list("alpha0"=alpha0, "alpha12"=alpha12.wmat)


# -- Only important in case variables are generated on which the network can be estimated
# mu: qxp matrix of weighting for mean
mu=matrix(0,q,p)
sweight=seq(-2.5,2.5,0.5)
mu[,sample(1:p, round(p*0.6))] <- sample(sweight, round(p*0.6)*q, replace = T)


# -- Save all relevant parameters in list
sparams <- list(iter=iter,n=n, p=p, q=q, alpha=alpha, mu=mu, 
                beta0=beta0, xbeta=xbeta, gbeta=gbeta, sthresh=dg.thresh, 
                eps.y=eps.y, eps.g=eps.g)



