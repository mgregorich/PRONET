#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Cockpit
#===============================================================================#



# ======================== Set up ==============================================
rm(list=ls())

# Functions
source("src/functions_main.R")
source("src/functions_aux.R")
source("src/setup.R")

# Paths
out.path <- "output/"
sim.path <- paste0(out.path, "sim_", Sys.Date(),"/")

if(!dir.exists(out.path)){dir.create(out.path)}
if(!dir.exists(sim.path)){dir.create(sim.path)}

# Save setup
setup <- readLines("src/setup.R")
write.table(setup, file = here::here(sim.path, "info_setup.txt"))

# ======================= Simulation ===========================================

# Barabasi-Albert model with quadratic preferential attachment for Bernoulli graph
BA.graph <- sample_pa(n=p, power=2, m=20, directed = F)                         # increase m to increase density

k=1
for(k in 1:nrow(scenarios)){
  print(paste0("Run scenario ",k,"/",nrow(scenarios)))
  scn <- scenarios[k,]
  filename <- paste0("sim_n",scn$n,"_p",scn$p,
                     "_beta",unlist(scn$beta.params)[1], "", unlist(scn$beta.params)[2],
                     "_DGt",scn$dg.thresh,"_epsY",scn$eps.y,"_epsG",scn$eps.g)
  
  simulate_pronet(iter = scn$iter,
                  n = scn$n, 
                  p = scn$p, 
                  q = scn$q, 
                  b0 = scn$b0, 
                  b1 = scn$b1, 
                  da.thresh = scn$da.thresh, 
                  dg.thresh = scn$dg.thresh,
                  beta.params = scn$beta.params,
                  alpha0.params = scn$alpha0.params,
                  alpha12.params = scn$alpha12.params,
                  X.params = scn$X.params,
                  eps.y = scn$eps.y, 
                  eps.g = scn$eps.g,
                  BA.graph = BA.graph,
                  filename = filename,
                  excel = scn$excel)
  if(scn$report){report_results(sim.path=sim.path, scen.nr=k, filename=filename)}
}



