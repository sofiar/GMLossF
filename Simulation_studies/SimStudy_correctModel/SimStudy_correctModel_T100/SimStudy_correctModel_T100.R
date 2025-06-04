###############################################################################
#
#
# Gompertz Model with Poisson Sampling Error: Correct Model
#
#
###############################################################################





############################### Initialization ################################

### set the seed
set.seed(1994)



##### set environment for simulation study

### compile stan model
# translation
rstan::stanc("../stan_model.stan")
# compilation and loading
stan_model = rstan::stan_model("../stan_model.stan")

### number of iterations in simulation study
niter = 500

### parameters for MCMC and mcmcEM
nsim = 1e+4
nsim_em = 1e+3
max_nsim_em = 2e+4
maxit_em = 1e+3

### sample size
T = 100



##### get true parameter values
theta_case1 = c(2, 0.22, -0.22)
theta_case2 = c(2, 0.22, -0.50)





################## Simulation Study #####################

##### parallel environment setup

### get cluster
cluster = parallel::makeCluster(24)

### get true number of iterations according to cluster size
niter = ceiling(niter / length(cluster)) * length(cluster)

### export variables
parallel::clusterExport(cl = cluster, varlist = ls(all.names = TRUE))

### identify cores
for (id_cluster in 1:length(cluster)) {
  parallel::clusterExport(cl = cluster[id_cluster], "id_cluster")
}

### create folder for reports file
suppressWarnings(dir.create("./reports"))

### open report file
parallel::clusterEvalQ(
  cl = cluster, sink(paste("./reports/report", id_cluster, ".txt", sep = ""))
)

### set seed
parallel::clusterSetRNGStream(cl = cluster)



##### do iterations in parallel
source("../SimStudy_correctModel_parexe.R")



##### parallel environment shut down

### close report file
parallel::clusterEvalQ(cl = cluster, sink())

### leave cluster
parallel::stopCluster(cl = cluster)



##### merge parallel executions
source("../SimStudy_correctModel_mergepar.R")





################################### Closure ###################################
save(
  list = ls(), file = paste("SimStudy_correctModel_T", T, ".RData", sep = "")
)
