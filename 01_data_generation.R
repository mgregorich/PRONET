############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data generation (code from https://github.com/shanghongxie/Covariate-adjusted-network)
############################################


generate_data <- function (n, p, q, zeta, alpha, obeta0, delta,beta0,xbeta, gbeta) {
  # n=params$n; p=params$p; q=params$q; 
  # zeta=params$zeta; alpha=params$alpha; delta=params$delta;
  # obeta0=params$obeta0; beta0=params$beta0;xbeta=params$xbeta; gbeta = params$gbeta
  
  
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
  
  df <- data.frame("Y"=Y, "X"=X, "FMI"=t(FMI), "MI"=MI)
  
  return(df)   
}


# tmp <- GenDataD(n=500, p=25, q=2, zeta, alpha, obeta0, delta,beta0,xbeta,mbeta,oeta)
# plot(tmp$FMI, tmp$Y)


