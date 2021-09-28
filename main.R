############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study
############################################

# 1. Generate data by conditional graphical model
# 2. estimate model by gaussian model


rm(list=ls())
library(mvtnorm)      # for multivariate normally distributed variables
library(igraph)       # Barabasi and Albert model
library(Rcpp)
library(RcppEigen)


##############################
###  Algorithm parameters  ###
##############################

n=500; 

## p: number of biomarker nodes; q: number of covariates; po: number of undirected edges; delta: variance of covariate xi
p=5; q=5; po=(p-1)*p/2;delta=1


### precision matrix
Omega=matrix(0,q,po)
Omega[1:3,1]=-c(0.5,1,1.5)
Omega[1:3,5]=-c(1,0.5,1.5)
Omega[1:3,8]=-c(-1.5,0.5,1)
Omega[1:3,10]=-c(0.5,1.5,-1)

obeta=as.vector(Omega)

## kappa: mean 

Kappa=matrix(0,q,p)
Kappa[1:3,1]=c(0.5,1,1.5)
Kappa[1:3,2]=c(1,0.5,1.5)
Kappa[1:3,3]=c(-1.5,1,0.5)

kbeta=as.vector(Kappa)

### 1/sigma^2
obeta0=rep(5,p)   ### sigma0^2=0.1

### 2nd stage parameter
beta0=0   ### intercept
xbeta=c(1,2,rep(0,q-2))  ## coefficients for covariate X
mbeta=c(1,3,rep(0,p-2))  ## coefficients for biomaker nodes M
oeta=c(1,0,0,0,2,0,0,0,0,0,rep(0,po-10))  ## coefficients for connections


GenDataD <- function (n, p, q, kbeta, obeta, obeta0, delta,beta0,xbeta,mbeta,oeta) {
  
  ## number of undirected edges
  po=(p-1)*p/2
  
  Omega=matrix(obeta,nrow=q)  # covariance between covariates
  Kappa=matrix(kbeta,nrow=q)  # mean of covariates
  
  ###  Generate X, M, MI  ###
  i=1; X=NULL; M=NULL; MI=NULL;idv=0;dv=0;imax=0;iboth=0
  
  repeat {
    dv=dv+1
    
    xi=mvrnorm(1,mu=rep(0,q),diag(delta,q))
    
    ox=c(xi%*%Omega); kx=c(xi%*%Kappa)
    
    Omegai=diag(obeta0); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        Omegai[s,r]=-ox[sr]; Omegai[r,s]=-ox[sr]
        sr=sr+1
      }
    }
    
    
    temmi=min(eigen(Omegai)$values)
    if (max(Omegai)>obeta0[1]) {imax=imax+1}
    if (max(Omegai)>obeta0[1] & (temmi<0.05)) {iboth=iboth+1}
    if (temmi<0.05) {idv=idv+1; next}
    
    Sigmai=solve(Omegai)
    
    ### mean matrix
    mui=Sigmai%*%kx
    
    
    ### generate biomarker nodes M
    mi=mvrnorm(1, mui, Sigmai)
    
    ### generate 2nd stage data ###
    ## Mutual Information
    count=0
    mii=numeric((p-1)*p/2); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        pho=ox[sr]^2/obeta0[s]/obeta0[r]
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
  
  
  # ###  Scaled for 2nd stage ###
  index=which(apply(MI,2,sum)!=0)
  inM=MI[,index]
  
  sdMI=matrix(0,n,po)
  temM=t(inM)-apply(inM,2,mean)
  temM=temM/sqrt(apply(temM^2,1,mean))
  temM=t(temM)
  sdMI[,index]=temM
  
  tem=t(X)-apply(X,2,mean)
  tem=tem/sqrt(apply(tem^2,1,mean))
  tem=t(tem)
  sdX=tem
  
  tem=t(M)-apply(M,2,mean)
  tem=tem/sqrt(apply(tem^2,1,mean))
  tem=t(tem)
  sdM=tem
  
  ####  Generate Y  ###
  Y=NULL
  for (i in 1:n) {
    xb=beta0+sum(sdX[i,]*xbeta)+sum(sdM[i,]*mbeta)+sum(sdMI[i,]*oeta)
    Y[i]=rnorm(1,mean=xb,sd=1)
  }
  
  
  return(list(X=X, M=M, Y=Y, MI=MI, sdX=sdX, sdM=sdM, sdMI=sdMI, idv=idv,dv=dv,count=count,imax=imax,iboth=iboth))   
}

M=tmp$M; X=tmp$X; Y=tmp$Y; MI=tmp$MI
sdM=tmp$sdM; sdX=tmp$sdX; sdMI=tmp$sdMI

