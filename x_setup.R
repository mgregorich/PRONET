# ============================================================================ #
# Author: MG
# Date: 07.02.2021
# Info: Simulation study - Setup
# ============================================================================ #


# ------------- Packages -------------------------
pacman::p_load(mvtnorm, igraph, NetworkToolbox, Rcpp, RcppEigen, MASS, lqmm, 
               stringr, future.apply, parallel, dplyr, tidyr, knitr, reshape2,
               refund, broom, cvTools, concreg)


# ------------- Code ----------------------------
out.path <- "../Output/"
sim.path <- paste0(out.path, "sim_", Sys.Date(),"/")
if(!dir.exists(sim.path)){dir.create(sim.path)}


## ----------- Parameters -------------------------
set.seed(666)

iter=5
n=250
q=2; 
delta=1                                                                    # q: number of covariates; delta: variance of covariate xi; qstar: number of latent processes
p=50;  po=(p-1)*p/2                                                             # p: number of biomarker nodes;  po: number of undirected edges
sthresh = 0.25                                                                   # Sparsification threshold for data gen
thresh.seq = seq(0,1,0.05)                                                     # Sparsification sequence for data ana

## mu: qxp matrix of weighting for mean 
mu=matrix(0,q,p)
sweight=seq(-2.5,2.5,0.5)
mu[,sample(1:p, round(p*0.6))] <- sample(sweight, round(p*0.6)*q, replace = T)

### omega: weighting matrix for the covariates to define the precision matrix
BA.graph <- sample_pa(n=p, power=2, m=15, directed = F)
BA.strc <- as.matrix(as_adjacency_matrix(BA.graph))

distr.params=list("beta"=c("shape1"=5, "shape2"=2), "unif"=c("min"=0, "max"=1), 
                  "alpha0.norm"=c("mean"=1, "sd"=0.5), "X.norm"=c("mean"=0, "sd"=1)) 
eta.params <- calc_eta_mean_and_var(norm.pars=distr.params$alpha0.norm, unif.pars=distr.params$unif)

alpha0 <- BA.strc[lower.tri(BA.strc)]
alpha0[alpha0==1] <- rnorm(sum(alpha0), distr.params$alpha0.norm[1], distr.params$alpha0.norm[2])
omega.imat=matrix(alpha0,1,po, byrow = T)

alpha12 <- rep(BA.strc[lower.tri(BA.strc)],q)
alpha12[alpha12==1] <- runif(sum(alpha12), distr.params$unif[1], distr.params$unif[2])                             
alpha12.wmat=matrix(alpha12,q,po, byrow = T)
alpha=list("alpha0"=alpha0, "alpha12"=alpha12.wmat)

## number of possible undirected edges
po=(p-1)*p/2

### 2nd stage parameter
beta0=1               # intercept
xbeta=2.5              # coefficients for covariate X
gbeta=3              # coefficients for network features fmi


sparams <- list(iter=iter,n=n, p=p, q=q, alpha=alpha, delta=delta, mu=mu, beta0=beta0, xbeta=xbeta, gbeta=gbeta, thresh=sthresh)

