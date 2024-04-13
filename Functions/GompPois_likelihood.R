################################################################################
#
# Gompertz Model with Poisson Sampling Error: Likelihood Computation
#
################################################################################


############### description
# This function computes the likelihood of the following model:
#
#           Nstar_t | N_t ~ Poisson(N_t)
#           N_t = exp(Z_t)
#           Z_t = a + (1 + b) Z_{t - 1} + eps_t
#           eps_t ~iid Gauss(0, sigmaSq)
#
# Thus, the likelihood function is pi(Nstar | a, b, sigmaSq). The loglikelihood
# is computed through an importance sampling algorithm. The importance density
# is constructed via a sequence of Gaussian distributions minimizing the
# Kolmogorov-Smirnov distance. Notice that a simple close formula exists for
# the parameters of the Gaussian.
#
# It is used an alternative parametrization, say:
#
#     theta1 = - a / b
#     theta2 = - sigmaSq / (b * (2 + b))
#
# which implies
#
#     Z_1, Z_2, ..., Z_T ~ Gauss_T ( theta1 1_T, theta2 SigmaBar)
#     SigmaBar_ij = (1 + b ) ^ abs(i - j)

#source("extraFunctions.R")

############### function
GompPois_likelihood = function(
  Nstar, theta1, theta2, b,
  nsim = 1e+5, log = FALSE, details = FALSE
) {

  # Nstar is the vector of the data, say Nstar_t, t = 1,2,...,T;
  # theta1 whose likelihood must be computed;
  # theta2 whose likelihood must be computed; 
  # b whose likelihood must be computed;
  # nsim is the number of simulations for the importance sampling;
  # log is a logical variable. Should the log-likelihood be computed?
  # details is a logical variable. If TRUE then also the weights and draws are
  #  reported.


  ###############--------------- function script 

  ########## initialization

  ### sample size
  T = length(Nstar)

  ### get values of alternative parametrization
  # location
  a = - theta1 * b
  # scale
  sigmaSq = - theta2 * b * (2 + b)

  ### matrix of drawn values
  Z = matrix(0, nrow = T, ncol = nsim)

  ### matrix of the conditional prior mean
  mean_Z = matrix(0, nrow = T, ncol = nsim)

  ### matrix of the conditional mean for importance density (IS)
  mean_Zis = matrix(0, nrow = T, ncol = nsim)

  ### matrix of the conditional standard for importance density (IS)
  sd_Zis = matrix(0, nrow = T, ncol = nsim)

  ### matrix of the log-weights
  logWeights = double(nsim)




  ########## importance sampling

  ##### start with t = 1

  ### prior mean
  mean_Z[1, ] = theta1

  ### get forward looking variable
  if(Nstar[2] != 0) frwrd_look = log(Nstar[2]) else frwrd_look = -1

  ### pre-mean for IS (mu)
  mean_Zis[1, ] = (
    theta1 * sigmaSq + (1 + b) * (frwrd_look - a) * theta2
  ) / (
    sigmaSq + (1 + b)^2 * theta2
  )

  ### pre-sd for IS (tauSq) and sd for IS (omegaSq)
  sd_Zis[1, ] = sqrt(sigmaSq * theta2 / (sigmaSq + (1 + b)^2 * theta2))

  ### mean for IS (xi)
  mean_Zis[1, ] = Nstar[1] * sd_Zis[1, ]^2 + mean_Zis[1, ] - lambertW_expArg(
    2 * log(sd_Zis[1, ]) + Nstar[1] * sd_Zis[1, ]^2 + mean_Zis[1, ]
  )

  ### sample from importance density
  Z[1, ] = mean_Zis[1, ] + sd_Zis[1, ] * rnorm(nsim)

  ### compute logweights
  logWeights = dpois(Nstar[1], lambda = exp(Z[1, ]), log = TRUE) +
    dnorm(Z[1,], mean = mean_Z[1, ], sd = sqrt(theta2), log = TRUE) -
    dnorm(Z[1,], mean = mean_Zis[1, ], sd = sd_Zis[1, ], log = TRUE)



  ##### continue up to T - 1
  for (t in 2:(T - 1)) {

    ### prior mean
    mean_Z[t, ] = a + (1 + b) * Z[t - 1, ]

    ### get forward looking variable
    if(Nstar[t + 1] != 0) frwrd_look = log(Nstar[t + 1]) else frwrd_look = -1

    ### pre-mean for IS (mu)
    mean_Zis[t, ] = (
      a + (1 + b) * (frwrd_look + Z[t - 1, ] - a)
    ) / (
      1 + (1 + b)^2
    )

    ### pre-sd for IS (tauSq) and sd for IS (omegaSq)
    sd_Zis[t, ] = sqrt(sigmaSq / (1 + (1 + b)^2))

    ### mean for IS (xi)
    mean_Zis[t, ] = Nstar[t] * sd_Zis[t, ]^2 + mean_Zis[t, ] - lambertW_expArg(
      2 * log(sd_Zis[t, ]) + Nstar[t] * sd_Zis[t, ]^2 + mean_Zis[t, ]
    )

    ### sample from importance density
    Z[t, ] = mean_Zis[t, ] + sd_Zis[t, ] * rnorm(nsim)

    ### compute logweights
    logWeights = logWeights +
      dpois(Nstar[t], lambda = exp(Z[t, ]), log = TRUE) +
      dnorm(Z[t, ], mean = mean_Z[t, ], sd = sqrt(sigmaSq), log = TRUE) -
      dnorm(Z[t, ], mean = mean_Zis[t, ], sd = sd_Zis[t, ], log = TRUE)
  }



  ##### end with t = T

  ### prior mean
  mean_Z[T, ] = a + (1 + b) * Z[T - 1, ]

  ### pre-mean for IS (mu)
  mean_Zis[T, ] = mean_Z[T, ]

  ### pre-sd for IS (tauSq) and sd for IS (omegaSq)
  sd_Zis[T, ] = sqrt(sigmaSq)

  ### mean for IS (xi)
  mean_Zis[T, ] = Nstar[T] * sd_Zis[T, ]^2 + mean_Zis[T, ] - lambertW_expArg(
    2 * log(sd_Zis[T, ]) + Nstar[T] * sd_Zis[T, ]^2 + mean_Zis[T, ]
  )

  ### sample from importance density
  Z[T, ] = mean_Zis[T, ] + sd_Zis[T, ] * rnorm(nsim)

  ### compute logweights
  logWeights = logWeights +
  dpois(Nstar[T], lambda = exp(Z[T, ]), log = TRUE) +
    dnorm(Z[T,], mean = mean_Z[T, ], sd = sqrt(sigmaSq), log = TRUE) -
    dnorm(Z[T,], mean = mean_Zis[T, ], sd = sd_Zis[T, ], log = TRUE)




  ########## finalization

  ### take the maximum log-weight
  maxLogWeight = max(logWeights)

  ### compute exponential of scaled log-weights
  expSLogWeights = exp(logWeights - maxLogWeight)

  ### compute logarithm of the unbiased estimator
  res = maxLogWeight - log(nsim) + log(sum(expSLogWeights))

  ### get the exponential of res if log is FALSE
  if (!log) res = exp(res)

  ### get results
  # detailed version
  if (details) {
    return(list(
      res = res,
      normWeights = expSLogWeights / sum(expSLogWeights),
      Z = Z
    ))
    # otherwise only the point value of the function
  } else {
    return(res)
  }

}





############### appendix

if(FALSE){ # appendix script is not run

  # debugging
  {
    rm(list = ls()); gc(); cat("\14")
    source("GompPois_likelihood.R")
    set.seed(1994)
    nsim = 1e+5
    T = 100
    theta1 = 1.9244
    theta2 = 0.4726 ^ 2 # that is 0.2233508
    b = -0.24
    log = TRUE
    details = TRUE

    Nstar = GompPois_simulate(T, theta1 = theta1, theta2 = theta2, b = b)
    start_time = Sys.time()
    res = GompPois_likelihood(
      Nstar, theta1 = theta1, theta2 = theta2, b = b, nsim = nsim,
      log = TRUE, details = details
    )
    end_time = Sys.time()
  }

  difftime(end_time, start_time)
  res$res
  plot(res$normWeights)
  1 / sum(res$normWeights^2)


  ###############
}