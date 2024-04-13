###############################################################
##################### Simulation analysis #####################
###############################################################
rm(list = ls(all.names = TRUE)); gc(); cat("\14")


source_dir = "Functions"

files = list.files(source_dir, pattern = "\\.R$", full.names = TRUE)
for (ifun in files) source(ifun)

set.seed(99)

# set paramters
theta1_s = 1.9244
theta2_s = 0.4726 ^ 2 # that is 0.2233508
b_s = c(-0.8, -0.24)
T = 100

a_s = -b_s * theta1_s
sigmaSq_s = -b_s * (2 + b_s) * theta2_s

niters = 2500
nsim = 1e+4
psi1 = 0.01
psi2 = 0.01
eta1 = 0
eta2 = 100
c = 1
starter = NULL
burn = 1
thin = 1
verbose = +Inf

probs = c(1:99) * 0.01

for (b_curr in b_s) {
  for (theta2_curr in theta2_s) {
    for (theta1_curr in theta1_s) {

      keep_percentiles = array(NA, dim = c(niters, length(probs), 3))
      keep_postMeans = matrix(NA, nrow = niters, ncol = 3)
      keep_time = double(niters)

      for (iiters in 1:niters) {
        #theta1_curr = theta1_s[1]
        #theta2_curr = theta2_s[1]
        #b_curr = b_s[1]; verbose = 100; iiters = 1
        Nstar = GompPois_simulate(
          T, theta1 = theta1_curr, theta2 = theta2_curr, b = b_curr
        )
        start_time = Sys.time()
        res = GompPois_quasiGibbs(
          nsim, Nstar = Nstar, psi1 = psi1, psi2 = psi2, eta1 = eta1,
          eta2 = eta2, c = c, starter = starter, burn = burn, thin = thin,
          verbose = verbose
        )
        keep_time[iiters] = difftime(Sys.time(), start_time, units = "sec")
        #print(difftime(Sys.time(), start_time, units = "sec"))
        keep_percentiles[iiters, , 1] = quantile(res$theta1, probs = probs)
        keep_percentiles[iiters, , 2] = quantile(res$theta2, probs = probs)
        keep_percentiles[iiters, , 3] = quantile(res$b, probs = probs)
        keep_postMeans[iiters, 1] = mean(res$theta1)
        keep_postMeans[iiters, 2] = mean(res$theta2)
        keep_postMeans[iiters, 3] = mean(res$b)
      }

      save(
        list = c("keep_time", "keep_percentiles", "keep_postMeans"),
        file = paste(
          "./Rdata/","postComp_", "theta1_", theta1_curr, "___", "theta2_",
          theta2_curr, "___", "b_", abs(b_curr), ".RData", sep = ""
        )
      )
    }
  }
}

# if (FALSE) {
#   rm(list = ls())
#   load("postComp_theta1_1.9244___theta2_0.22335076___b_0.8.RData")
#   keep_percentiles[, 1, 1]
#   keep_postMeans[, 1]
#   keep_time
# }