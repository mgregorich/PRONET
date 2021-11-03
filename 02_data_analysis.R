############################################
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Data analysis 
############################################

rm(list=ls())

# -------------  Data generation
source("01_data_generation.R")


# -------------  Data analysis
# Network data
data.network <- data.sim[,str_detect(colnames(data.sim), "MI.")]

# CC for threshold sequence
list.gvars <- lapply(1:nrow(data.network), function(x) evalSSN(eweights=data.network[x,], msize=p))
data.gvars <- data.frame(t(do.call(cbind, list.gvars)))
colnames(data.gvars) <- paste0("T",seq(0,1,0.05))

# Add outcome Y
data.gvars <- cbind("Y"=data.sim$Y, data.gvars)

# Pick model with best RMSE
res <- data.frame(t(sapply(colnames(data.gvars)[-1], function(x) evalLM(x, dt=data.gvars[,c("Y", x)]))))
res
which.min(res$rmse)

