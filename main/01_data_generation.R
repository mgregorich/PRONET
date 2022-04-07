#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data generation (code from https://github.com/shanghongxie/Covariate-adjusted-network)
#===============================================================================#


generate_data <- function (n, p, q, mu, alpha, distr.params, eta.params, obeta0, delta, beta0, xbeta, gbeta, eps, sthresh) {
  # n=sparams$n; p=sparams$p; q=sparams$q; mu=sparams$mu; distr.params=distr.params;
  # alpha=sparams$alpha; delta=sparams$delta;obeta0=sparams$obeta0; beta0=sparams$beta0;xbeta=sparams$xbeta;
  # gbeta = sparams$gbeta; eps = list("y"=sparams$eps.y, "g"=sparams$eps.g);
  # sthresh=sparams$sthresh
  
  # -------  Network generation
  data.graph <- genIndivNetwork(n=n, p=p, q=q, alpha=alpha, 
                                distr.params=distr.params, eta.params = eta.params,
                                delta, mu=mu)
  
  # -------  Compute network features
  GE <- abs(data.graph$GE)
  GE.thres <- t(sapply(1:nrow(GE), function(x){
    if(length(sthresh)==1){
      thr.weight=sthresh
    }else{
      thr.weight=sample(sthresh,1, replace=T) 
      }
    mat <- VecToSymMatrix(0, side.entries = GE[x,], mat.size = p)
    res.mat <- data.frame(weightThresholding(mat, w=thr.weight, method = "trim")$adj)
    res <- res.mat[upper.tri(res.mat)]
    }))
  GE.fea <- data.frame(t(apply(GE.thres, 1, function(x) calcGraphFeatures(VecToSymMatrix(0, x, p)))))

  
  # ------- Outcome generation
  Y=NULL
  for (i in 1:n) {
    xb=beta0+sum(GE.fea[i,1]*gbeta)
    Y[i]=rnorm(1,mean=xb,sd=eps$y)
  }
  GE.noisy = GE + matrix(rnorm(length(GE),0, sd=eps$g), nrow = 250)
  true.R2 = cor(Y, GE.fea[,1])^2
  
  df <- data.frame("Y"=Y, "X"=data.graph$X, "GE.fea"=GE.fea, "GE"=GE, "GE.noisy"=GE.noisy, "true.tau"=sthresh, "true.R2"=true.R2)
  
  return(df)   
}


