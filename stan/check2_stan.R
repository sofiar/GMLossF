setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rstan::stanc("stan_model.stan")
stan_model = rstan::stan_model("stan_model.stan")
# rstan::stan(
#   file = "stan_model.stan", data = list(
#     M = length(gse::AmericanRedstart2[1:5, "redstart"]),
#     Nstar = gse::AmericanRedstart2[1:5, "redstart"]
#   ),
#   chain = 1, cores = 1, iter = 10, warmup = 5
# )
# res_our = gse::GompPois_flatNIG_mcmc(
#   nsim = 1e+4, Nstar = gse::AmericanRedstart2[, "redstart"],
#   phi1 = 0.1, phi2 = 0.1, eta1 = 0, eta2 = 100,
# )
rstan::stan

set.seed(1994)

y = gse::GompPois_rng(T = 30, theta1 = 1.90, theta2 = 0.2, b = -0.2)


res_stan = rstan::sampling(
  object = stan_model, data = list(M = length(y), Nstar = y),
  chain = 1, cores = 5, iter = 1e+4 + 1e+3, warmup = 1e+3
)
time_stan = rstan::get_elapsed_time(res_stan)
# 1560.38 seconds (Total), T = 100, cores = 1
# 1213.97 seconds (Total), T = 100, cores = 5

sample_stan = rstan::extract(res_stan, permuted = FALSE)

summ_stan = rstan::summary(res_stan)

summ_stan$summary[1:3, "n_eff"]
sns::ess(sample_stan[, 1, c("theta1", "theta2", "b")])
sns::ess(sample_stan[, 1, c("theta1", "theta2", "b")]) / sum(time_stan)
sns::ess(sample_stan[, 1, c("theta1", "theta2", "b")]) / 1560.38
sns::ess(sample_stan[, 1, c("theta1", "theta2", "b")]) / 1213.97
plot(sample_stan[, 1, "b"], type = "l")

time_our = Sys.time()
res_our = gse::GompPois_flatNIG_mcmc(
  nsim = 1e+4, Nstar = y, phi1 = 0.1, phi2 = 0.1, eta1 = 0, eta2 = 100,
)
time_our = difftime(Sys.time(), time_our, units = "secs")
# 32 secs, T = 100
sns::ess(t(res_our)[, c("theta1", "theta2", "b")])
sns::ess(t(res_our)[, c("theta1", "theta2", "b")]) / as.numeric(time_our)
sns::ess(t(res_our)[, c("theta1", "theta2", "b")]) / 32

plot(res_our["theta2", ], res_our["b", ], type = "p", cex = 0.5, xlim = c(0, 5))

save(list = ls(), file = "check2_stan.RData")

colMeans(sample_stan[, 1, c("theta1", "theta2", "b")])