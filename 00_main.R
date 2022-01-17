#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Cockpit
#===============================================================================#


rm(list=ls())



# ---------- Set up --------------------------------

pacman::p_load(mvtnorm, igraph, NetworkToolbox, Rcpp, RcppEigen, MASS, lqmm, 
               stringr, future.apply, parallel, dplyr, tidyr, knitr, reshape2)

source("functions_aux.R")
source("01_data_generation.R", print.eval=TRUE)
source("02_data_analysis.R", print.eval=TRUE)
source("03_summarize_results.R", print.eval=TRUE)

out.path <- "../Output/"
sim.path <- paste0(out.path, "sim_", Sys.Date(),"/")
dir.create(sim.path)


## ----------- Parameters -------------------------
set.seed(1234)

iter=100
n=100
q=2; delta=1                                                                    # q: number of covariates; delta: variance of covariate xi; qstar: number of latent processes
p=25;  po=(p-1)*p/2                                                             # p: number of biomarker nodes;  po: number of undirected edges
sthresh = 0.25                                                                   # Sparsification threshold for data gen
thresh.seq = seq(0,1,0.025)                                                     # Sparsification sequence for data ana

## mu: qxp matrix of weighting for mean 
mu=matrix(0,q,p)
sweight=seq(-2.5,2.5,0.5)
mu[,sample(1:p, round(p*0.6))] <- sample(sweight, round(p*0.6)*q, replace = T)

### omega: weighting matrix for the covariates to define the precision matrix
omega.distr=list("icpt"=c("par1"=1.15, "par2"=0.75), "weights"=c("par1"=0, "par2"=1)) 
omega.strc <- as.matrix(as_adjacency_matrix(sample_pa(n=25, power=1, m=4)))
omega.icpt <- omega.strc[lower.tri(omega.strc)]
omega.icpt[omega.icpt==1] <- rbeta(sum(omega.icpt), omega.distr$icpt[1], omega.distr$icpt[2])*8
omega.imat=matrix(omega.icpt,1,po, byrow = T)

omega.weights <- rep(omega.strc[lower.tri(omega.strc)],q)
omega.weights[omega.weights==1] <- rnorm(sum(omega.weights), omega.distr$weights[1], omega.distr$weights[2])                              
omega.wmat=matrix(omega.weights,q,po, byrow = T)
omega=list("icpt"=omega.icpt, "wmat"=omega.wmat)

## number of possible undirected edges
po=(p-1)*p/2

### 2nd stage parameter
beta0=10               # intercept
xbeta=2.5              # coefficients for covariate X
gbeta=2.5              # coefficients for network features fmi


sparams <- list(iter=iter,n=n, p=p, q=q, omega=omega, delta=delta, mu=mu, beta0=beta0, xbeta=xbeta, gbeta=gbeta, thresh=sthresh)


# ------------ SIMULATION -------------------

wrapper_sim <- function(sparams, tseq){
  
  data.iter <- generate_data(n=sparams$n, p=sparams$p, q=sparams$q, 
                             omega=sparams$omega, delta=sparams$delta,mu=sparams$mu,
                             obeta0=sparams$obeta0, beta0=sparams$beta0,xbeta=sparams$xbeta, gbeta = sparams$gbeta)
  results.iter <- analyse_data(data.iter, tseq=tseq)
  
  return(results.iter)
}

plan(multisession, workers=detectCores()/2)
results.sim <- future_lapply(1:iter, function(x) wrapper_sim(sparams, tseq = thresh.seq), future.seed = 666)
results.sim <- do.call(rbind, results.sim)
plan(sequential)

output.sim <- summarize_results(res=results.sim, sparams, tseq=thresh.seq)
output.sim


# -- One run for testing purposes
# data.test <- generate_data(n=sparams$n, p=sparams$p, q=sparams$q,
#                            mu=sparams$mu, omega=sparams$omega, delta=sparams$delta,
#                            beta0=sparams$beta0,xbeta=sparams$xbeta, gbeta = sparams$gbeta)
# results.test <- analyse_data(data.test, tseq=thresh.seq)


# ---------- Markdown Report ----------------
# rmarkdown::render(
#   "04_report_results.Rmd",
#   params = list(sparams = sparams, tbl_list=output.sim),
#   output_file = paste0(sim.path,"Report_" ,Sys.Date(), ".html")
# )


