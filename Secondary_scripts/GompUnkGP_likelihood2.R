##########################################################################################
#
# Gompertz Model with Unknown Response II Version: Gaussian White
#  Noise and Poisson Sampling Distribution
# 
##########################################################################################


############### description
# This function computes the likelihood of the following model:
#
#           Nstar_t | N_t ~ Poisson(N_t)
#           N_t = log(Z_t)
#           Z_t = a + (1+b) Z_{t-1} + eps_t
#           eps_t ~iid Gauss(0, sigmaSq)
#
# Thus, the likelihood function is pi(Nstar | a, b, sigmaSq). It is computed through 
# importance sampling.
#
# It is used an alternative parametrization, say:
#
#           theta1 = - a / b
#           theta2 = - sigmaSq / ( b (2+b) )
#
# which implies
#
#           Z_1, Z_2, ..., Z_T ~ Gauss_T ( theta1 1_T, theta2 SigmaBar)
#           SigmaBar_ij = (1+b) ^ abs(i-j)



############### function
GompUnkGP_likelihood2 = function(
    theta1, theta2, b, Nstar,
    nsim, log = FALSE, details = FALSE
){
  
  # theta1 whose likelihood must be computed;
  # theta2 whose likelihood must be computed; 
  # b whose likelihood must be computed;
  # Nstar is the vector of the observed variables, say Nstar_t, t = 1,2,...,T;
  # nsim is the number of simulations for the importance sampling;
  # log is a logical variable. Should the log-likelihood be computed?
  # details is a logical variable. If TRUE then also the weights and draws are reported.
  
  
  ###############--------------- function script 
  
  ########## initialization
  
  ### sample size
  T = length(Nstar)
  
  ### get values of alternative parametrization
  # location
  a = - theta1 * b
  # scale
  sigmaSq = - theta2 * b*(2+b)
  
  ### matrix of drawn values
  Z = matrix(nrow = T, ncol = nsim)
  
  ### matrix of the conditional prior mean
  mean_Z = matrix(nrow = T, ncol = nsim)
  
  ### matrix of the conditional mean (LA)
  mean_Zla = matrix(nrow = T, ncol = nsim)
  
  ### matrix of the conditional standard deviation (LA)
  sd_Zla = matrix(nrow = T, ncol = nsim)
  
  ### matrix of the log-weights
  logWeights = vector(length = nsim) 
  
  ### adjusted values of logarithm of observed data
  # take the logarithm
  logNstar_adj = log(Nstar)
  # adjust
  logNstar_adj[ is.infinite(logNstar_adj) ] = 0.006737947
  
  
  
  
  ########## importance sampling
  
  ##### start with t = 1
  
  ### prior mean
  mean_Z[1,] = theta1
  
  ### get parameters of LA
  # mean
  mean_Zla[1,] = (
    theta1*sigmaSq + (1+b) * (logNstar_adj[2]-a) * theta2
  ) / (
    sigmaSq + (1+b)^2 * theta2
  )
  # standard deviation
  sd_Zla[1,] = sqrt(
    (
      sigmaSq*theta2
    )/(
      sigmaSq + (1+b)^2 * theta2
    )
  )
  
  ### draw Z_1
  Z[1,] = rnorm( nsim, mean = mean_Zla[1,], sd = sd_Zla[1,] )
  
  ### compute the log-weights
  logWeights = dpois(Nstar[1], lambda = exp(Z[1,]), log = TRUE) +
    dnorm(Z[1,], mean = theta1, sd = sqrt(theta2), log = TRUE) -
    dnorm(Z[1,], mean = mean_Zla[1,], sd = sd_Zla[1,], log = TRUE)
  
  
  
  ##### continue up to T-1
  for( t in 2:(T-1) ){
    
    ### prior mean
    mean_Z[t,] = a + (b+1) * Z[t-1,]
    
    ### get parameters of LA
    # mean
    mean_Zla[t,] = (
      a + (1+b) * (logNstar_adj[t+1]+Z[t-1,]-a)
    ) / (
      1 + (1+b)^2
    )
    # standard deviation
    sd_Zla[t,] = sqrt(
      (
       sigmaSq
      ) / (
        sigmaSq + (1+b)^2
      )
    )
    
    ### draw Z_t-1
    Z[t,] = rnorm( nsim, mean = mean_Zla[t,], sd = sd_Zla[t,] )
    
    ### compute the log-weights
    logWeights = logWeights + dpois(Nstar[t], lambda = exp(Z[t,]), log = TRUE) +
      dnorm(Z[t,], mean = mean_Z[t,], sd = sqrt(sigmaSq), log = TRUE) -
      dnorm(Z[t,], mean = mean_Zla[t,], sd = sd_Zla[t,], log = TRUE)
    
  }
  
  
  
  ##### end with T
  ### prior mean
  mean_Z[T,] = a + (b+1) * Z[T-1,]
  
  ### get parameters of LA
  # mean
  mean_Zla[T,] = a + (b+1) * Z[T-1,]
  # standard deviation
  sd_Zla[T,] = sqrt(sigmaSq)
  
  ### draw Z_T
  Z[T,] = rnorm( nsim, mean = mean_Zla[T,], sd = sd_Zla[T,] )
  
  ### compute the log-weights
  logWeights = logWeights + dpois(Nstar[T], lambda = exp(Z[T,]), log = TRUE) +
    dnorm(Z[T,], mean = mean_Z[T,], sd = sqrt(sigmaSq), log = TRUE) -
    dnorm(Z[T,], mean = mean_Zla[T,], sd = sd_Zla[T,], log = TRUE)
  
  
  
  
  ########## finalization
  
  ### take the maximum log-weight
  maxLogWeight = max(logWeights)
  
  ### compute logarithm of the unbiased estimator
  res = maxLogWeight - log(nsim) + log( sum(
    exp( logWeights - maxLogWeight )
  ) )
  
  ### get the exponential of res if log is FALSE
  if ( !log ) { res = exp(res) }
  
  ### get results
  # detailed version
  if(details){
    return(list(
      res = res,
      logWeights = logWeights,
      Z = Z
    ))
    # otherwise only the point value of the function
  }else{
    return(res)
  }
  
}





############### appendix

if(FALSE){ # appendix script is not run
  
  # debugging
  rm(list = ls(all.names = TRUE)); gc(); cat("\14")
  set.seed(1994)
  theta1 = 0
  theta2 = 1
  b = -0.8
  Nstar = rpois(n = 500, lambda = exp(theta1))
  nsim = 1e+4
  log = TRUE
  details = TRUE
  trial = GompUnkGP_likelihood2(theta1, theta2, b, Nstar, nsim, log, details)
  
  plot(trial$logWeights, cex = 0.5, pch = 19)
  plot(
    exp(trial$logWeights-trial$res) / sum(exp(trial$logWeights-trial$res)),
    cex = 0.5, pch = 19
  )
  
  
  ###############
  
  
  # trials  
  
}