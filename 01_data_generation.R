#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data generation (code from https://github.com/shanghongxie/Covariate-adjusted-network)
#===============================================================================#


generate_data <- function (n, p, q, mu, alpha, distr.params, eta.params, obeta0, delta,beta0, xbeta, gbeta, eps=0.1) {
  # n=sparams$n; p=sparams$p; q=sparams$q; mu=sparams$mu; distr.params=distr.params;
  # alpha=sparams$alpha; delta=sparams$delta;obeta0=sparams$obeta0; beta0=sparams$beta0;xbeta=sparams$xbeta; gbeta = sparams$gbeta
  
  # -------  Network generation
  data.graph <- genIndivNetwork(n=n, p=p, q=q, alpha=alpha, 
                                distr.params=distr.params, eta.params = eta.params,
                                delta, mu=mu)
  
  # -------  Compute network features
  GE <- abs(data.graph$GE)
  GE.thresh <- data.frame(weightThresholding(GE, w=0.25, method = "trim")$adj)
  GE.fea <- data.frame(t(apply(GE.thresh, 1, function(x) calcGraphFeatures(VecToSymMatrix(0, x, p)))))

  #--------- Outcome generation
  Y=NULL
  for (i in 1:n) {
    xb=beta0+sum(GE.fea[i,1]*gbeta)
    Y[i]=rnorm(1,mean=xb,sd=eps)
  }
  
  df <- data.frame("Y"=Y, "X"=data.graph$X, "GE.fea"=GE.fea, "GE"=GE)
  
  return(df)   
}


