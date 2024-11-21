###############################################################
##################### Simulation analysis #####################
#### True likelihood comparison: Bayesian vs frequentist   ####
###############################################################

rm(list = ls(all.names = TRUE)); gc(); cat("\14")

source_dir = '/u/ruizsuar/GMLossF/Functions'

files = list.files(source_dir, pattern = "\\.R$", full.names = TRUE)
for (ifun in files) source(ifun)

set.seed(1606)

# set parameters
theta1_s = 1.9244
theta2_s = 0.4726 ^ 2   
b_s = c(-0.24,-0.8) #c(-0.8, -0.24)
T_s = c(30,100,500)#c(30, 100, 500)  # 30  100 - 500

a_s = -b_s * theta1_s
sigmaSq_s = -b_s * (2 + b_s) * theta2_s

niters = 600 #600
nsim = 1e+4

phi1 = 0.5
phi2 = 0.5
zeta1 = 0
zeta2 = 1
nu = 2
c = 1
starter = NULL
burn = 1
thin = 1
verbose = +Inf

probs = c(1:99) * 0.01

nscenarios = length(b_s) * length(theta2_s) * length(theta1_s) * length(T_s)
iscenarios = 0

print(paste('starting time is ',Sys.time(),sep=''))
      
for (T_curr in T_s) {
  for (b_curr in b_s) {
    for (theta2_curr in theta2_s) {
      for (theta1_curr in theta1_s) {

        keep_percentiles_bay = array(NA, dim = c(niters, length(probs), 3))
        keep_postMeans_bay = matrix(NA, nrow = niters, ncol = 3)
        keep_time_bay = double(niters)

        keep_percentiles_eb = array(NA, dim = c(niters, length(probs), 3))
        keep_postMeans_eb = matrix(NA, nrow = niters, ncol = 3)
        keep_time_eb = double(niters)
        
        keep_mle_freq = matrix(NA, nrow = niters, ncol = 3)
        keep_time_freq = double(niters)
        
        iscenarios = iscenarios + 1

        for (iiters in 1:niters) {
          
          ###############################################################
          ######################## Simulate data ########################
          ###############################################################
          
          Nstar = GompPois_simulate(
            T = T_curr, theta1 = theta1_curr, theta2 = theta2_curr, b = b_curr
          )
          

          ###############################################################
          ################### Fit full Bayesian model ###################
          ###############################################################
          
          start_time = Sys.time()
          # res_bay = GompPois_UScaledPriors_quasiGibbs (
          # nsim, Nstar = Nstar, phi1 = phi1, phi2 = phi2, nu = nu,
          # zeta1 = zeta1, zeta2 = zeta2, c = c, starter = starter,
          # burn = burn, thin = thin, verbose = verbose
          # )
          res_bay = GompPois_ScaledPriors_quasiGibbs (
          nsim, Nstar = Nstar, phi1 = phi1, phi2 = phi2, nu = nu,
          zeta1 = zeta1, zeta2 = zeta2, c = c, starter = starter,
          burn = burn, thin = thin, verbose = verbose
          )

          # save output
          keep_time_bay[iiters] = difftime(Sys.time(), start_time, units = "sec")
          #print(difftime(Sys.time(), start_time, units = "sec"))
          keep_percentiles_bay[iiters, , 1] = quantile(res_bay$theta1, probs = probs)
          keep_percentiles_bay[iiters, , 2] = quantile(res_bay$theta2, probs = probs)
          keep_percentiles_bay[iiters, , 3] = quantile(res_bay$b, probs = probs)
          keep_postMeans_bay[iiters, 1] = mean(res_bay$theta1)
          keep_postMeans_bay[iiters, 2] = mean(res_bay$theta2)
          keep_postMeans_bay[iiters, 3] = mean(res_bay$b)
        
          ####################################################################
          ################### Fit Empirical Bayes model ######################
          ####################################################################

          start_time = Sys.time()
          outputEB = GompPois_hyperTuner(Nstar = Nstar, maxit = 1e+3, 
          abstol = sqrt(.Machine$double.eps),nsim = 1e+4, starter = NULL,
          verbose = verbose)

          phi1_est = outputEB[3]
          phi2_est = outputEB[4]
          zeta1_est = outputEB[1]
          zeta2_est = outputEB[2]
          
          res_eb = GompPois_quasiGibbs (
          nsim, Nstar = Nstar, psi1 = phi1_est, psi2 = phi2_est,
          eta1 = zeta1_est, eta2 = zeta2_est, c = c, starter = starter,
          burn = burn, thin = thin, verbose = verbose
          ) 

          # save output
          keep_time_eb[iiters] = difftime(Sys.time(), start_time, units = "sec")
          #print(difftime(Sys.time(), start_time, units = "sec"))
          keep_percentiles_eb[iiters, , 1] = quantile(res_eb$theta1, probs = probs)
          keep_percentiles_eb[iiters, , 2] = quantile(res_eb$theta2, probs = probs)
          keep_percentiles_eb[iiters, , 3] = quantile(res_eb$b, probs = probs)
          keep_postMeans_eb[iiters, 1] = mean(res_eb$theta1)
          keep_postMeans_eb[iiters, 2] = mean(res_eb$theta2)
          keep_postMeans_eb[iiters, 3] = mean(res_eb$b)

          ####################################################################
          #################### Fit Freq approach ######################
          ####################################################################
          
          start_time = Sys.time()
          res_freq = GompPois_mcmc_em(Nstar = Nstar, maxit = 1e+3, nsim = 1e+4, 
          verbose =verbose)
          
          keep_time_freq[iiters] = difftime(Sys.time(), start_time, units = "sec")
          keep_mle_freq[iiters,1] = res_freq[1] # theta1
          keep_mle_freq[iiters,2] = res_freq[2] # theta2
          keep_mle_freq[iiters,3] = res_freq[3] # b

          ####################################################################

          print(paste('iteration ',iiters,' of ', niters,' for scenario ', iscenarios,
          ' of ', nscenarios, ' at time ', Sys.time(),sep=''))
          
          # Save RData 
          save(
            list = c("keep_time_bay", "keep_percentiles_bay", "keep_postMeans_bay",
                    "keep_time_eb", "keep_percentiles_eb", "keep_postMeans_eb",
                    "keep_mle_freq","keep_time_freq","niters","nsim","phi1",
                    "phi2","zeta2","nu","c","starter","burn","thin","verbose"),
            file = paste(
              "/u/ruizsuar/GMLossF/Rdata/","TLs_BayVsFreq", "theta1_", theta1_curr, "___", "theta2_",
              theta2_curr, "___", "b_", abs(b_curr), "___","T_", T_curr, ".RData", sep = "" 
            )
          )
        
        }
     
    }
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

# save.image(file='/u/ruizsuar/GMLossF/Rdata/Error_Environment_Oct16.RData')