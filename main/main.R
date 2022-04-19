#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Cockpit
#===============================================================================#



# ======================== Set up ==============================================
rm(list=ls())

# Functions
source("main/functions_main.R")
source("main/functions_aux.R")
source("main/setup.R")


# Paths
out.path <- "output/"
sim.path <- paste0(out.path, "sim_", Sys.Date(),"/")
if(!dir.exists(out.path)){dir.create(out.path)}
if(!dir.exists(sim.path)){dir.create(sim.path)}


# ======================= Simulation ===========================================
# varying parameters: beta shape parameter, dg threshold, eps_y

simulate_pronet(main.params)
report_results(sim.path)



