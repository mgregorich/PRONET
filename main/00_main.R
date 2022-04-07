#===============================================================================#
# Author: MG
# Date: 06.07.2021
# Info: Simulation study - Cockpit
#===============================================================================#


rm(list=ls())


# ======================== Set up ==============================================
source("main/x_functions.R")
source("main/x_setup.R")
source("main/01_data_generation.R", print.eval=TRUE)
source("main/02_data_analysis.R", print.eval=TRUE)

out.path <- "output/"
sim.path <- paste0(out.path, "sim_", Sys.Date(),"/")
if(!dir.exists(out.path)){dir.create(out.path)}
if(!dir.exists(sim.path)){dir.create(sim.path)}

# ======================= Simulation ===========================================
# -- Data generation & analysis
plan(multisession, workers=detectCores()/2)
results.sim <- results.iter <-list()
results.sim <- future_lapply(1:iter, function(x){
  data.iter <- generate_data(n=sparams$n, p=sparams$p, q=sparams$q,
                             alpha=sparams$alpha, delta=sparams$delta,mu=sparams$mu,
                             distr.params = distr.params, eta.params = eta.params,
                             obeta0=sparams$obeta0, beta0=sparams$beta0,xbeta=sparams$xbeta,
                             gbeta = sparams$gbeta,  eps = list("y"=sparams$eps.y, "g"=sparams$eps.g),
                             sthresh=sparams$sthresh)
  results.iter <- analyse_data(data.iter, true.params=true.params, tseq=da.thresh)
  return(results.iter)
  }, future.seed = 666)
plan(sequential)


# -- Summarize results
source("main/03_summarize_results.R")


# -- Report results
rmarkdown::render(
  "main/04_report_results.Rmd",
  params = list(output_dir=sim.path, sim_dir=paste0("../",sim.path)),
  output_dir = output_dir,
  output_file = paste0("Report_" ,Sys.Date(), ".html"))
browseURL(file.path(paste0(sim.path,"Report_" ,Sys.Date(), ".html")))




# --------------- One run for testing purposes -----------------
# results.test <- list()
# data.test <- generate_data(n=sparams$n, p=sparams$p, q=sparams$q,
#                            mu=sparams$mu, alpha=sparams$alpha, delta=sparams$delta,
#                            distr.params=distr.params, eta.params = eta.params,
#                            beta0=sparams$beta0,xbeta=sparams$xbeta, gbeta = sparams$gbeta,
#                            eps = list("y"=sparams$eps.y, "g"=sparams$eps.g),
#                            sthresh=sparams$sthresh)
# results.test <- analyse_data(data.test, tseq=da.thresh, true.tau=sparams$sthresh, k=5)

