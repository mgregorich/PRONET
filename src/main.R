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
if(dir.exists(sim.path)){invisible(do.call(file.remove, list(list.files(sim.path, full.names = T))))
}else{dir.create(sim.path)}

# Save setup
setup <- readLines("src/setup.R")
write.table(setup, file = here::here(sim.path, "info_setup.txt"))


# ======================= Simulation ===========================================

k=1
for(k in 1:nrow(scenarios)){
  print(paste0("Run scenario ",k,"/",nrow(scenarios)))
  scn <- scenarios[k,]
  file_dgthresh = names(scn$dg.thresh)
  filename <- paste0("sim_i", scn$iter,"_n",scn$n,"_p",scn$p,
                     "_beta",unlist(scn$beta.params)[1], "", unlist(scn$beta.params)[2],
                     "_DGt",file_dgthresh,"_epsY",scn$eps.y,"_epsG",scn$eps.g)
  
  simulate_pronet(iter = scn$iter,
                  n = scn$n, 
                  p = scn$p, 
                  q = scn$q, 
                  b0 = scn$b0, 
                  b1 = scn$b1, 
                  dg.thresh = scn$dg.thresh,
                  beta.params = scn$beta.params,
                  alpha0.params = scn$alpha0.params,
                  alpha12.params = scn$alpha12.params,
                  X1.params = scn$X1.params,
                  X2.params = scn$X2.params,
                  eps.y = scn$eps.y, 
                  eps.g = scn$eps.g,
                  BA.graph = BA.graph,
                  filename = filename,
                  excel = scn$excel)
}


# Summarize all scenarios
sim.files <- list.files(sim.path, pattern = "sim_")
sim.all <- lapply(sim.files, function(x) cbind_results(x, sim.path))
tbl_scens <- do.call(rbind, sim.all)

saveRDS(tbl_scens, here::here(sim.path, "tbl_scenarios_results.rds"))  
write.xlsx(tbl_scens, here::here(sim.path, "tbl_scenarios_results.xlsx"), overwrite = T)

# Generate Markdown report with results
# report_simresults(sim.path, filename=paste0("report_results_", Sys.Date()))
