###############################################################
##################### Simulation analysis #####################
## Comparison True likelihood with Composite Likelihood one ###
###############################################################

rm(list = ls(all.names = TRUE)); gc(); cat("\14")

source_dir = '/u/ruizsuar/GMLossF/Functions'

files = list.files(source_dir, pattern = "\\.R$", full.names = TRUE)
for (ifun in files) source(ifun)

#set.seed(99)
set.seed(1606)

# set paramters
theta1_s = 1.9244
theta2_s = 0.4726 ^ 2 # that is 0.2233508
b_s = c(-0.8, -0.24)
T = 30

a_s = -b_s * theta1_s
sigmaSq_s = -b_s * (2 + b_s) * theta2_s

niters = 10
nsim = 1e+4
phi1 = 0.5
phi2 = 0.5
zeta1 = 0
zeta2 = 1
nu = 2
c = 2
starter = NULL
burn = 1
thin = 1
verbose = +Inf

probs = c(1:99) * 0.01

ntriplets = length(b_s) * length(theta2_s) * length(theta1_s)
itriplet = 0

print(paste('starting time is ',Sys.time(),sep=''))
      

for (b_curr in b_s) {
  for (theta2_curr in theta2_s) {
    for (theta1_curr in theta1_s) {

      keep_percentiles_tl = array(NA, dim = c(niters, length(probs), 3))
      keep_postMeans_tl = matrix(NA, nrow = niters, ncol = 3)
      keep_time_tl = double(niters)

      keep_percentiles_cl = array(NA, dim = c(niters, length(probs), 3))
      keep_postMeans_cl = matrix(NA, nrow = niters, ncol = 3)
      keep_time_cl = double(niters)
      keep_acceptance_rate = double(niters)

      itriplet = itriplet + 1

      for (iiters in 1:niters) {
        
        ###############################################################
        ######################## Simulate data ########################
        ###############################################################

        Nstar = GompPois_simulate(
          T, theta1 = theta1_curr, theta2 = theta2_curr, b = b_curr
        )
        start_time = Sys.time()
        
        ###############################################################
        ################ Fit true likelihood model ####################
        ###############################################################
        
        res_tl = GompPois_UScaledPriors_quasiGibbs (
         nsim, Nstar = Nstar, phi1 = phi1, phi2 = phi2, nu = nu,
         zeta1 = zeta1, zeta2 = zeta2, c = c, starter = starter,
         burn = burn, thin = thin, verbose = verbose
        ) 

        # save output
        keep_time_tl[iiters] = difftime(Sys.time(), start_time, units = "sec")
        #print(difftime(Sys.time(), start_time, units = "sec"))
        keep_percentiles_tl[iiters, , 1] = quantile(res_tl$theta1, probs = probs)
        keep_percentiles_tl[iiters, , 2] = quantile(res_tl$theta2, probs = probs)
        keep_percentiles_tl[iiters, , 3] = quantile(res_tl$b, probs = probs)
        keep_postMeans_tl[iiters, 1] = mean(res_tl$theta1)
        keep_postMeans_tl[iiters, 2] = mean(res_tl$theta2)
        keep_postMeans_tl[iiters, 3] = mean(res_tl$b)
      

        ###############################################################
        ############## Fit composite likelihood model #################
        ###############################################################

        # get Omega 
        post_sample_likelihood = matrix(nrow = 3, ncol = nsim)
        post_sample_likelihood[1, ] = res_tl$theta1
        post_sample_likelihood[2, ] = res_tl$theta2
        post_sample_likelihood[3, ] = res_tl$b

        ### to avoid infinite values
        post_sample_likelihood[2, post_sample_likelihood[2, ] == 0] = 0 + .Machine$double.eps
        post_sample_likelihood[3, post_sample_likelihood[3, ] == -1] = -1 + .Machine$double.eps
        post_sample_likelihood[3, post_sample_likelihood[3, ] == 0] = 0 - .Machine$double.eps

        transf_post_sample = post_sample_likelihood
        transf_post_sample[2 ,] = log(post_sample_likelihood[2, ])
        transf_post_sample[3, ] = log(
         -post_sample_likelihood[3, ] / (1 + post_sample_likelihood[3, ])
         )
        delta = diff(t(transf_post_sample))
        Omega = t(delta) %*% delta / (nsim - 1)

        # Fit model 
        res_cl = GompPois_UScaledPriors_compositeLike(
        nsim, Nstar = Nstar, phi1 = phi1, phi2 = phi2, nu = nu,
        zeta1 = zeta1, zeta2 = zeta2, c = c, Omega = Omega, starter = starter, 
        burn = burn, thin = thin, verbose = verbose
       ) 

      # save the output        
        keep_percentiles_cl[iiters, , 1] = quantile(res_cl$theta1, probs = probs)
        keep_percentiles_cl[iiters, , 2] = quantile(res_cl$theta2, probs = probs)
        keep_percentiles_cl[iiters, , 3] = quantile(res_cl$b, probs = probs)
        keep_postMeans_cl[iiters, 1] = mean(res_cl$theta1)
        keep_postMeans_cl[iiters, 2] = mean(res_cl$theta2)
        keep_postMeans_cl[iiters, 3] = mean(res_cl$b)
        keep_acceptance_rate[iiters] = res_cl$true_accep_rate

      print(paste('iteration ',iiters,' of ', niters,' for triplet ',itriplet,
      ' of triplet ',ntriplets, ' at time ',Sys.time(),sep=''))
      
      
      }
     # Save RData 
      save(
        list = c("keep_time_tl", "keep_percentiles_tl", "keep_postMeans_tl",
                 "keep_time_cl", "keep_percentiles_cl", "keep_postMeans_cl",
                 "keep_acceptance_rate","niters","nsim","phi1","phi2","zeta2",
                 "nu","c","starter","burn","thin","verbose"),
        file = paste(
          "/u/ruizsuar/GMLossF/Rdata/","SimStudy_UscaledPrior", "theta1_", theta1_curr, "___", "theta2_",
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