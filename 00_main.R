#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Cockpit
#===============================================================================#


rm(list=ls())


# ---------- Set up --------------------------------
source("x_functions.R")
source("x_setup.R")

# ------------ SIMULATION -------------------
source("01_data_generation.R", print.eval=TRUE)
source("02_data_analysis.R", print.eval=TRUE)

plan(multisession, workers=detectCores()/2)
results.sim <- results.iter <-list()
results.sim <- future_lapply(1:iter, function(x){
  data.iter <- generate_data(n=sparams$n, p=sparams$p, q=sparams$q,
                             alpha=sparams$alpha, delta=sparams$delta,mu=sparams$mu,
                             distr.params = distr.params, eta.params = eta.params,
                             obeta0=sparams$obeta0, beta0=sparams$beta0,xbeta=sparams$xbeta, 
                             gbeta = sparams$gbeta,  eps = sparams$eps, sthresh=sparams$sthresh)
  results.iter <- analyse_data(data.iter, tseq=thresh.seq)
  return(results.iter)
  }, future.seed = 666)
plan(sequential)


# ---------- Summarise results --------------
source("03_summarize_results.R")

# ---------- Markdown Report ----------------
rmarkdown::render(
  "04_report_results.Rmd",
  output_file = paste0(sim.path,"Report_" ,Sys.Date(), ".html")
)




# # --------------- One run for testing purposes -----------------
# results.test <- list()
# data.test <- generate_data(n=sparams$n, p=sparams$p, q=sparams$q,
#                            mu=sparams$mu, alpha=sparams$alpha, delta=sparams$delta,
#                            distr.params=distr.params, eta.params = eta.params,
#                            beta0=sparams$beta0,xbeta=sparams$xbeta, gbeta = sparams$gbeta, eps = sparams$eps,
#                            sthresh=sparams$sthresh)
# results.test <- analyse_data(data.test, tseq=thresh.seq, k=5)

