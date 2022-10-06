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

# --- Run through all scenarios
# plan(multisession, workers = detectCores()*.75)
# invisible(future_lapply(1:nrow(scenarios), function(k) simulate_scenario(scn=scenarios[k,]), future.seed = TRUE))
# plan(sequential)

# # --- Run through parts
# scenarios.partial <- scenarios[1:360,]
# plan(multisession, workers = detectCores()*.75)
# future_lapply(1:nrow(scenarios.partial), function(k) run_scenario(scn=scenarios.partial[k,]), future.seed = TRUE)
# plan(sequential)
# 
scenarios.partial <- scenarios[361:405,]
plan(multisession, workers = detectCores()*.5)
future_lapply(1:nrow(scenarios.partial), function(k) run_scenario(scn=scenarios.partial[k,]), future.seed = TRUE)
plan(sequential)


# --- Summarize all scenarios
sim.files <- list.files(sim.path, pattern = "sim_")
sim.all <- evaluate_scenarios(sim.files, sim.path)
saveRDS(sim.all, here::here(sim.path, "tbl_scenario_results.rds"))  


# --- Generate Markdown report with results
# report_simresults(sim.path, filename=paste0("report_results_2022-08-09"))
