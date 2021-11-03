############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data generation (code from https://github.com/shanghongxie/Covariate-adjusted-network)
############################################


rm(list=ls())

pacman::p_load(mvtnorm, igraph, NetworkToolbox, Rcpp, RcppEigen, MASS, lqmm, stringr, future.apply)
source("functions_aux.R")
set.seed(666)


## ----------- Parameters -------------------------

# n: sample size
n=500; 
# q: number of covariates; delta: variance of covariate xi; qstar: number of latent processes
q=2; delta=1
# p: number of biomarker nodes;  po: number of undirected edges
p=25;  po=(p-1)*p/2
# Sparsification threshold
sthres = 0.25

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
beta0=10      ### intercept
xbeta=2.5    ## coefficients for covariate X
gbeta=2.5   ## coefficients for network features fmi


## ------------ Data gen alg -------------------------
GenDataD <- function (n, p, q, zeta, alpha, obeta0, delta,beta0,xbeta,mbeta,oeta) {
  
  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, M, MI for each subject i separately ###
  i=1; X=NULL; M=NULL; MI=NULL; FMI=NULL
  
  # ---- (1) Network Generation
  repeat {

    # Covariates X_i
    xi=MASS::mvrnorm(1,mu=rep(0,q),diag(delta,q))
    
    # Omega_i: Precision matrix, Kappa_i:mean
    ox=c(xi%*%alpha); kx=c(xi%*%zeta)
    Omegai=VecToSymMatrix(obeta0, -ox)
    
    # No covariance matrix is generated that is singular
    if(!is.positive.definite(Omegai)){Omegai<-make.positive.definite(Omegai, tol=1e-3)}
    
    # Covariance matrix: Inverse of Omegai
    Sigmai=solve(Omegai)
    
    ### mean matrix ???
    #mui=Sigmai%*%kx
    
    ### generate biomarker nodes M
    #mi=MASS::mvrnorm(1, mui, Sigmai)
    mi=MASS::mvrnorm(1, kx, Sigmai)
    
    ## Partial correlation - Network edge weights
    count=0; sr=1
    mii=numeric((p-1)*p/2); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        # rho^2=pho
        pho=-ox[sr]/sqrt(obeta0[s]*obeta0[r])
        if (abs(pho)>0.9999){count=count+1}
        pho=ifelse(abs(pho)>0.9999, sign(pho)*0.9999, pho)
        mii[sr]=pho
        sr=sr+1
      }
    }

    X=rbind(X,xi)      # covariates
    M=rbind(M,mi)      # network nodes
    MI=rbind(MI,mii)   # network edges
    
    if (i==n) break
    i=i+1
  }
  
  # ------ (2) Network features
  #### Sparsification ###
  MI = abs(MI)
  MI.thres = (MI > sthres)*1
  FMI <- data.frame(t(apply(MI.thres,1, function(x) calcGraphFeatures(VecToSymMatrix(0, x, p)))))

  #------ (3) Outcome generation
  Y=NULL
  for (i in 1:n) {
    #xb=beta0+sum(X[i,]*xbeta)+sum(M[i,]*mbeta)+sum(MI[i,]*oeta)
    xb=beta0+sum(FMI[i]*gbeta)
    Y[i]=rnorm(1,mean=xb,sd=0.25)
  }
  
  return(list(X=X, Y=Y, M=M, MI=MI, FMI=FMI))   
}

tmp <- GenDataD(n, p, q, zeta, alpha, obeta0, delta,beta0,xbeta,mbeta,oeta)


data.sim <- data.frame("Y"=tmp$Y, "X"=tmp$X, "FMI"=t(tmp$FMI), "MI"=tmp$MI)

plot(data.sim$FMI, data.sim$Y)


