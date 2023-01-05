#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Cockpit
#===============================================================================#



# ======================== Set up ==============================================
rm(list=ls())

# Load & save setup
source("src/setup.R")

if(!dir.exists(here::here("output"))){dir.create(here::here("output"))}
if(dir.exists(sim.path)){invisible(do.call(file.remove, list(list.files(sim.path, full.names = T))))
}else{dir.create(sim.path)}

setup <- readLines("src/setup.R")
write.table(setup, file = here::here(sim.path, "info_setup.txt"))



# ======================= Simulation ===========================================

# --- Run through all scenarios
plan(multisession, workers = detectCores()*.75)
invisible(future_lapply(1:nrow(scenarios), function(k) simulate_scenario(scn=scenarios[k,]), future.seed=0xBEEF))
plan(sequential)



# --- Summarize all scenarios
sim.files <- list.files(sim.path, pattern = "sim_")
sim.all <- evaluate_scenarios(sim.files, sim.path)
saveRDS(sim.all, here::here(sim.path, "tbl_scenario_results.rds"))  


# --- Generate Markdown report with results
report_simresults(sim.path, filename=paste0("report_results_", sim.date))
