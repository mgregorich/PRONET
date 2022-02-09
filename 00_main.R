#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Cockpit
#===============================================================================#


rm(list=ls())


# ---------- Set up --------------------------------
source("x_setup.R")
source("x_functions.R")

# ------------ SIMULATION -------------------
source("01_data_generation.R", print.eval=TRUE)
source("02_data_analysis.R", print.eval=TRUE)

plan(multisession, workers=detectCores()/2)
results.sim <- future_lapply(1:iter, function(x) wrapper_sim(sparams, tseq = thresh.seq), future.seed = 666)
results.sim <- do.call(rbind, results.sim)
plan(sequential)

# # -- One run for testing purposes
# data.test <- generate_data(n=sparams$n, p=sparams$p, q=sparams$q,
#                            mu=sparams$mu, omega=sparams$omega, delta=sparams$delta,
#                            beta0=sparams$beta0,xbeta=sparams$xbeta, gbeta = sparams$gbeta)
# results.test <- analyse_data(data.test, tseq=thresh.seq)


# ---------- Summarise results --------------
source("03_summarize_results.R")


# ---------- Markdown Report ----------------
rmarkdown::render(
  "04_report_results.Rmd",
  output_file = paste0(sim.path,"Report_" ,Sys.Date(), ".html")
)


