############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data generation (code from https://github.com/shanghongxie/Covariate-adjusted-network)
############################################


rm(list=ls())

pacman::p_load(mvtnorm, igraph, Rcpp, RcppEigen, MASS)


## ----------- Parameters -------------------------

# n: sample size
n=500; 
# q: number of covariates; delta: variance of covariate xi
q=5; delta=1
# p: number of biomarker nodes;  po: number of undirected edges
p=5;  po=(p-1)*p/2

### alpha: weighting matrix for the covariates to define the precision matrix
alpha=matrix(0,q,po)
alpha[1:3,1]=-c(0.5,1,1.5)
alpha[1:3,5]=-c(1,0.5,1.5)
alpha[1:3,8]=-c(-1.5,0.5,1)
alpha[1:3,10]=-c(0.5,1.5,-1)

## zeta: qxp matrix of weighting for mean 
zeta=matrix(0,q,p)
zeta[1:3,1]=c(0.5,1,1.5)
zeta[1:3,2]=c(1,0.5,1.5)
zeta[1:3,3]=c(-1.5,1,0.5)

### 1/sigma^2
obeta0=rep(5,p)   ### sigma0^2=0.2

### 2nd stage parameter
beta0=0   ### intercept
xbeta=c(1,2,rep(0,q-2))  ## coefficients for covariate X
mbeta=c(1,3,rep(0,p-2))  ## coefficients for biomaker nodes M
oeta=c(1,0,0,0,2,0,0,0,0,0,rep(0,po-10))  ## coefficients for connections


## ------------ Data gen alg -------------------------
GenDataD <- function (n, p, q, zeta, alpha, obeta0, delta,beta0,xbeta,mbeta,oeta) {
  
  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, M, MI for each subject i separately ###
  i=1; X=NULL; M=NULL;  M2=NULL;MI=NULL;idv=0;dv=0;imax=0;iboth=0
  
  repeat {
    dv=dv+1
    
    # Covariates X_i
    xi=MASS::mvrnorm(1,mu=rep(0,q),diag(delta,q))
    
    # Omega_i: Precision matrix, Kappa_i:mean
    ox=c(xi%*%alpha); kx=c(xi%*%zeta)
   
    Omegai=diag(obeta0); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        Omegai[s,r]=-ox[sr]; Omegai[r,s]=-ox[sr]
        sr=sr+1
      }
    }
    
    #???
    temmi=min(eigen(Omegai)$values)
    if (max(Omegai)>obeta0[1]) {imax=imax+1}
    if (max(Omegai)>obeta0[1] & (temmi<0.05)) {iboth=iboth+1}
    if (temmi<0.05) {idv=idv+1; next} # ensures no covariance matrix is generated that is singular
    
    
    # Covariance matrix: Inverse of Omegai
    Sigmai=solve(Omegai)
    
    ### mean matrix ???
    #mui=Sigmai%*%kx
    
    ### generate biomarker nodes M
    #mi=MASS::mvrnorm(1, mui, Sigmai)
    mi=MASS::mvrnorm(1, kx, Sigmai)
    
    ### generate 2nd stage data ###
    ## Mutual Information
    count=0
    mii=numeric((p-1)*p/2); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        # rho^2=pho
        pho=ox[sr]^2/(obeta0[s]*obeta0[r])
        if (abs(pho)>0.9999){
          count=count+1
        }
        
        pho=ifelse(abs(pho)>0.9999, sign(pho)*0.9999, pho)
        mii[sr]=-log(1-pho)/2
        sr=sr+1
      }
    }
    
    X=rbind(X,xi)
    M=rbind(M,mi)
    MI=rbind(MI,mii)
    
    if (i==n) break
    i=i+1
  }
  
  #### Sparsification ###

  ####  Generate Y  ###
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

