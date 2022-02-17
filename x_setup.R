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
set.seed(1234)

iter=10
n=250
q=2; delta=1                                                                    # q: number of covariates; delta: variance of covariate xi; qstar: number of latent processes
p=50;  po=(p-1)*p/2                                                             # p: number of biomarker nodes;  po: number of undirected edges
sthresh = 0.25                                                                   # Sparsification threshold for data gen
thresh.seq = seq(0,1,0.025)                                                     # Sparsification sequence for data ana

## mu: qxp matrix of weighting for mean 
mu=matrix(0,q,p)
sweight=seq(-2.5,2.5,0.5)
mu[,sample(1:p, round(p*0.6))] <- sample(sweight, round(p*0.6)*q, replace = T)

### omega: weighting matrix for the covariates to define the precision matrix
omega.distr=list("icpt"=c("par1"=1.15, "par2"=0.2), "weights"=c("par1"=1, "par2"=0.5)) 
BA.graph <- sample_pa(n=p, power=2, m=15, directed = F)
BA.strc <- as.matrix(as_adjacency_matrix(BA.graph))
omega.icpt <- BA.strc[lower.tri(BA.strc)]
omega.icpt[omega.icpt==1] <- trunc_rnorm(sum(omega.icpt), omega.distr$icpt[1], omega.distr$icpt[2])
omega.imat=matrix(omega.icpt,1,po, byrow = T)

omega.weights <- rep(BA.strc[lower.tri(BA.strc)],q)
omega.weights[omega.weights==1] <- trunc_rnorm(sum(omega.weights), omega.distr$weights[1], omega.distr$weights[2], min = 0, max=2)                             
omega.wmat=matrix(omega.weights,q,po, byrow = T)
omega=list("icpt"=omega.icpt, "wmat"=omega.wmat)

## number of possible undirected edges
po=(p-1)*p/2

### 2nd stage parameter
beta0=1               # intercept
xbeta=2.5              # coefficients for covariate X
gbeta=3              # coefficients for network features fmi


sparams <- list(iter=iter,n=n, p=p, q=q, omega=omega, delta=delta, mu=mu, beta0=beta0, xbeta=xbeta, gbeta=gbeta, thresh=sthresh)

