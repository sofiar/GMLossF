library(rstan)

GompPois_rng = function(T, theta1, theta2, b) {
  
  ### initialize result object
  res = double(T)
  
  ### get the intercept
  a = -b * theta1
  
  ### get the conditional variance
  sigmaSq = -b * (2 + b) * theta2
  
  ### sample the first draw from the stationary distribution
  res[1] = theta1 + sqrt(theta2) * rnorm(1)
  
  ### sample the remaining draws from the conditional distribution
  for (t in 2:T) {
    res[t] = a + (1 + b) * res[t - 1] + sqrt(sigmaSq) * rnorm(1)
  }
  
  ### sample observed variables from Poisson
  res = rpois(T, lambda = exp(res))
  
  return(res)
  
}

stanc("stan_model.stan")
 
theta1_real = c(1.90, 1.90)
theta2_real = c(0.55, 0.2)
b_real = c(-0.2, -0.35)
nreps = 120

output_theta1_mean = matrix(NA,ncol=2, nrow=nreps)
colnames(output_theta1_mean) = c('triplet_Bay', 'triplet_freq')
output_theta2_mean = matrix(NA,ncol=2, nrow=nreps)
colnames(output_theta2_mean) = c('triplet_Bay', 'triplet_freq')
output_b_mean = matrix(NA,ncol=2, nrow=nreps)
colnames(output_b_mean) = c('triplet_Bay', 'triplet_freq')


for (s in 1:length(theta1_real)){
  
  for (rep in 1:nreps){
    
    Nstar = GompPois_rng (100, theta1 = theta1_real[s], theta2 = theta2_real[s], 
                          b = b_real[s]) 
    
    # Run Stan model
    ts_dat = list(M =length(Nstar), Nstar = Nstar)
    
    out = stan(file = 'stan_model.stan', data = ts_dat, chain=1, cores = 5,
                iter = 1e+04 )
    
    extracts = rstan::extract(out, permuted = TRUE) # return a list of arrays 
    
    output_theta1_mean[rep,s] = mean(extracts$theta1)
    output_theta2_mean[rep,s] = mean(extracts$theta2)
    output_b_mean[rep,s] = mean(extracts$b)
    
  }
}


	
save.image(file='StanNsamp100.RData')
