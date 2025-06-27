###############################################################################
#
#
# Gompertz Model with Poisson Sampling Error: American Redstart Two Dataset
#
#
###############################################################################





############################### Initialization ################################

### set the seed
set.seed(1994)



##### set environment for real data analysis

### get data
AmericanRedstart2 = gse::AmericanRedstart2

### compile stan model
# translation
rstan::stanc("./stan_model.stan")
# compilation and loading
stan_model = rstan::stan_model("./stan_model.stan")

### parameters for MCMC and mcmcEM
nsim = 1e+4
nsim_em = 1e+3
max_nsim_em = 2e+4
maxit_em = 1e+3




################## Real Data Analysis #####################

### maximum likelihood estimator
time = Sys.time()
res_mle = gse::GompPois_mcmcEM(
  Nstar = AmericanRedstart2$redstart, nsim = nsim_em, max_nsim = max_nsim_em,
  maxit = maxit_em, verbose = 1,
)
time_mle = difftime(Sys.time(), time, units = "secs")

### Bayesian approach (Gibbs)
time = Sys.time()
res_Gibbs = gse::GompPois_flatNIG_mcmc(
  nsim, Nstar = AmericanRedstart2$redstart,
  phi1 = 0.01, phi2 = 0.01, eta1 = 0, eta2 = 100, verbose = 1e+3
)
time_Gibbs = difftime(Sys.time(), time, units = "secs")

### Bayesian approach (stan)
object_stan = rstan::sampling(
  object = stan_model, data = list(
    M = length(AmericanRedstart2$redstart), Nstar = AmericanRedstart2$redstart
  ),
  chain = 1, cores = 1, iter = nsim + 1e+3, warmup = 1e+3,
  show_messages = FALSE, refresh = 0
)
time_stan = rstan::get_elapsed_time(object_stan)
sample_stan = rstan::extract(object_stan, permuted = FALSE)
res_stan = matrix(
  nrow = 3, ncol = nsim, dimnames = list(c("theta1", "theta2", "b"))
)
res_stan[1, ] = sample_stan[, 1, "theta1"]
res_stan[2, ] = sample_stan[, 1, "theta2"]
res_stan[3, ] = sample_stan[, 1, "b"]





################################### Closure ###################################
save(list = ls(), file = "RealData_AmericanRedstart2.RData")
