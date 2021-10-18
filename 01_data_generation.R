############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data generation (code from https://github.com/shanghongxie/Covariate-adjusted-network)
############################################


rm(list=ls())

pacman::p_load(mvtnorm, igraph, Rcpp, RcppEigen, MASS, lqmm)
source("functions_aux.R")

## ----------- Parameters -------------------------

# n: sample size
n=500; 
# q: number of covariates; delta: variance of covariate xi; qstar: number of latent processes
q=5; delta=1; qstar=1
# p: number of biomarker nodes;  po: number of undirected edges
p=5;  po=(p-1)*p/2

### alpha: weighting matrix for the covariates to define the precision matrix
alpha=matrix(0,q,po)
alpha[1:3,1]=c(0.5,1,1.5)
alpha[1:3,2]=c(1.5,1,-1)
alpha[1:3,5]=c(2.5,0.5,1.5)
alpha[1:3,8]=c(-3,0.5,1)
alpha[1:3,9]=c(-2,1.5,2)
alpha[1:3,10]=c(-1,-1.5,-1)
alpha = alpha[1:qstar,1:po]

## zeta: qxp matrix of weighting for mean 
zeta=matrix(0,q,p)
zeta[1:3,1]=c(0.5,1,0.5)
zeta[1:3,2]=c(1,0.5,1.5)
zeta[1:3,3]=c(-1.5,1,0.5)
zeta = zeta[1:qstar,1:p]

### 1/sigma^2
obeta0=rep(5,p)   ### sigma0^2=0.2

### 2nd stage parameter
beta0=0   ### intercept
xbeta=c(1,2,rep(0,q-qstar-2))  ## coefficients for covariate X
mbeta=c(1,3,rep(0,p-2))  ## coefficients for biomaker nodes M
oeta=c(1,0,0,0,2,0,0,0,0,0,rep(0,po-10))  ## coefficients for connections


## ------------ Data gen alg -------------------------
GenDataD <- function (n, p, q, zeta, alpha, obeta0, delta,beta0,xbeta,mbeta,oeta) {
  
  q=qstar
  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, M, MI for each subject i separately ###
  i=1; X=NULL; M=NULL;  M2=NULL;MI=NULL;idv=0;dv=0;imax=0;iboth=0
  
  # ---- (1) Network Generation
  repeat {
    dv=dv+1
    
    # Covariates X_i
    xi=MASS::mvrnorm(1,mu=rep(1,q),diag(delta,q))
    
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
    paste0(round(mii,3))
    
    X=rbind(X,xi)
    M=rbind(M,mi)
    MI=rbind(MI,mii)
    
    if (i==n) break
    i=i+1
  }
  
  # ------ (2) Network features
  #### Sparsification ###


  #------ (3) Outcome generation
  Y=NULL
  for (i in 1:n) {
    #xb=beta0+sum(X[i,]*xbeta)+sum(M[i,]*mbeta)+sum(MI[i,]*oeta)
    xb=beta0+sum(X[i,]*xbeta)+sum(MI[i,]*oeta)
    Y[i]=rnorm(1,mean=xb,sd=1)
  }
  
  return(list(X=X, M=M, Y=Y, MI=MI, idv=idv,dv=dv,count=count,imax=imax,iboth=iboth))   
}

tmp <- GenDataD(n, p, q, zeta, alpha, obeta0, delta,beta0,xbeta,mbeta,oeta)

M=tmp$M;
X=tmp$X; Y=tmp$Y; MI=tmp$MI
sdM=tmp$sdM; sdX=tmp$sdX; sdMI=tmp$sdMI

