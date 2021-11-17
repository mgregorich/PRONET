############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Cockpit
############################################


rm(list=ls())



# ---------- Set up --------------------------------

pacman::p_load(mvtnorm, igraph, NetworkToolbox, Rcpp, RcppEigen, MASS, lqmm, stringr, future.apply, parallel)

source("functions_aux.R")
source("01_data_generation.R", print.eval=TRUE)
source("02_data_analysis.R", print.eval=TRUE)
source("03_summarize_results.R", print.eval=TRUE)

out.path <- "../Output/"
sim.path <- paste0(out.path, "sim_", Sys.Date())
dir.create(sim.path)




## ----------- Parameters -------------------------
iter=5
n=500;                                                                          # n: sample size
q=2; delta=1                                                                    # q: number of covariates; delta: variance of covariate xi; qstar: number of latent processes
p=25;  po=(p-1)*p/2                                                             # p: number of biomarker nodes;  po: number of undirected edges
sthres = 0.25                                                                   # Sparsification threshold for data gen
thresh.seq = seq(0,1,0.05)                                                      # Sparsification sequence for data ana

### alpha: weighting matrix for the covariates to define the precision matrix
alpha=matrix(0,q,po)
sweight=seq(-2.5,2.5,0.5)
alpha[,sample(1:po, round(po*0.6))] <- sample(sweight,round(po*0.6)*q, replace = T)

## zeta: qxp matrix of weighting for mean 
zeta=matrix(0,q,p)
zeta[,sample(1:p, round(p*0.6))] <- sample(sweight, round(p*0.6)*q, replace = T)

### 1/sigma^2
obeta0=rep(10,p)   ### sigma0^2=0.1

### 2nd stage parameter
beta0=10               # intercept
xbeta=2.5              # coefficients for covariate X
gbeta=2.5              # coefficients for network features fmi

params <- list(iter=iter,n=n, p=p, q=q, zeta=zeta, alpha=alpha, obeta0=obeta0, delta=delta, beta0=beta0, xbeta=xbeta, gbeta=gbeta)




# ------------ SIMULATION -------------------

wrapper_sim <- function(params, tseq){
  
  data.iter <- generate_data(n=params$n, p=params$p, q=params$q, 
                             zeta=params$zeta, alpha=params$alpha, delta=params$delta,
                             obeta0=params$obeta0, beta0=params$beta0,xbeta=params$xbeta, gbeta = params$gbeta)
  results.iter <- analyse_data(data.iter, tseq=tseq)
  
  return(results.iter)
}

plan(multisession, workers=detectCores()/2)
results.sim <- future_lapply(1:iter, function(x) wrapper_sim(params, tseq = thresh.seq), future.seed = 666)
plan(sequential)

output.sim <- summarize_results(res=results.sim, params, tseq=thresh.seq)
output.sim


# data.sim <- generate_data(n=500, p=25, q=2, zeta, alpha, obeta0, delta,beta0,xbeta,mbeta,oeta)
# results <- analyse_data(data.sim)
