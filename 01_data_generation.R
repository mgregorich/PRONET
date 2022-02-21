#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data generation (code from https://github.com/shanghongxie/Covariate-adjusted-network)
#===============================================================================#


generate_data <- function (n, p, q, mu, omega, obeta0, delta,beta0,xbeta, gbeta) {
  # n=sparams$n; p=sparams$p; q=sparams$q; mu=sparams$mu; omega=sparams$omega; delta=sparams$delta;obeta0=sparams$obeta0; beta0=sparams$beta0;xbeta=sparams$xbeta; gbeta = sparams$gbeta
  
  
  ## number of possible undirected edges
  po=(p-1)*p/2
  
  ###  Generate X, GN, GE for each subject i separately ###
  i=1; X=NULL; GN=NULL; GE=NULL; GE.fea=NULL
  
  # ---- (1) Network Generation
  repeat {

    # Covariates X_i
    xi=MASS::mvrnorm(1,mu=rep(0,q),diag(delta,q))
    
    # Omega_i: Precision matrix, Kappa_i:mean
    omega.icpt = omega[[1]]
    omega.wmat = omega[[2]]
    ox = omega.icpt + c(xi%*%omega.wmat)
    obeta0 = rep(max(abs(ox)),p)   
    
    Mui=c(xi%*%mu)
    Omegai=VecToSymMatrix(obeta0, -ox)
    
    # No covariance matrix is generated that is singular
    if(!is.positive.definite(Omegai)){Omegai<-make.positive.definite(Omegai, tol=1e-3)}
    
    # Covariance matrix: Inverse of Omegai
    Sigmai=solve(Omegai)
    
    ### mean matrix ???
    #mui=Sigmai%*%Mui
    
    ### generate biomarker nodes M
    gn=MASS::mvrnorm(1, Mui, Sigmai)
    
    ## Partial correlation - Network edge weights
    sr=1
    ge=numeric((p-1)*p/2); sr=1
    for (s in 1:(p-1)) {
      for (r in (s+1):p) {
        pho=-ox[sr]/sqrt(obeta0[s]*obeta0[r])
        ge[sr]=pho
        sr=sr+1
      }
    }
    
    X=rbind(X,xi)      # covariates
    GN=rbind(GN,gn)    # graph nodes
    GE=rbind(GE,ge)   # graph edges
    
    if (i==n) break
    i=i+1
  }
  
  # ------ (2) Network features
  #### Sparsification ###
  GE = abs(GE)
  GE.thres = (GE > sparams$thresh)*1
  GE.fea <- data.frame(t(apply(GE.thres,1, function(x) calcGraphFeatures(VecToSymMatrix(0, x, p)))))

  #------ (3) Outcome generation
  Y=NULL
  for (i in 1:n) {
    #xb=beta0+sum(X[i,]*xbeta)+sum(M[i,]*mbeta)+sum(GE[i,]*oeta)
    xb=beta0+sum(GE.fea[i,1]*gbeta)
    Y[i]=rnorm(1,mean=xb,sd=0.1)
  }
  
  df <- data.frame("Y"=Y, "X"=X, "GE.fea"=GE.fea, "GE"=GE)
  
  return(df)   
}


# tmp <- GenDataD(n=500, p=25, q=2, mu, omega, obeta0, delta,beta0,xbeta,mbeta,oeta)
# plot(tmp$FGE, tmp$Y)


