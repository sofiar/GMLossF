################################################################################
#
# Gompertz Model with Poisson Sampling Error: Markov Chain Monte Carlo
#  Posterior Simulation Using Quasi Gibbs Sampler
#
################################################################################


############### description
# This function samples from the posterior distribution of the following model:
#
#           Nstar_t | N_t ~ Poisson(N_t)
#           N_t = exp(Z_t)
#           Z_t = a + (1 + b) Z_{t - 1} + eps_t
#           eps_t ~iid Gauss(0, sigmaSq)
#
# It is used an alternative parametrization, that is:
#
#     theta1 = - a / b
#     theta2 = - sigmaSq / (b * (2 + b))
#     beta = 1 / (1 + exp(b))
#
# which implies
#
#     Z_1, Z_2, ..., Z_T ~ Gauss_T ( theta1 1_T, theta2 SigmaBar)
#     SigmaBar_ij = (1 + b ) ^ abs(i - j)
#
# the priors are
#
#     beta ~ cauchy(0, c)
#     theta2 ~ invlomax(phi1, phi2)
#     theta1 | theta2 ~ student-t(nu, zeta1, zeta2 * theta2)
#
#
# The following scheme is used
#
#     (1) sample V, W, U | theta1, theta2, beta
#         (1.1) sample V | beta (posterior of a logistic Kolmogorov)
#         (1.2) sample U | theta2 (gamma)
#         (1.3) sample W | theta1, theta2 (inverse gamma)
#
#     (2) for t = 1, 2, ..., T
#         sample Z_t | Nstar_t, Z_{-t}, beta, theta2, theta1
#         (using acceptance-rejection)
#
#     (3) beta | Z, V, theta1, theta2 using elliptical slice sampling
#
#     (4) theta2, theta1 | beta, Z
#         (4.1) theta2 | beta, Z ~ invgamma(...)
#         (4.2) theta1 | theta2, beta, Z ~ gauss(...)

# source("extraFunctions.R")
# source("GompPois_likelihood.R")
# source("rpostlogiskolmo.R")

############### function
GompPois_UScaledPriors_quasiGibbs = function(
  nsim, Nstar, phi1, phi2, nu, zeta1, zeta2, c,
  starter = NULL, burn = 1, thin = 1, verbose = +Inf
) {

  # nsim is the number of simulations to be stored;
  # Nstar is the vector of the data, say Nstar_t, t = 1,2,...,T;
  # phi1, phi2, nu, zeta1, zeta2, c are the prior parameters, see introduction;
  # starter is a list providing the starter point;
  # burn is the the number of scans to be discarded as burn-in;
  # thin is the thinning parameter;
  # verbose is the period for printing status of the chain.

  # starter must be a list with the following structure:
  #  starter$theta1, starter$theta2, starter$b , starter$Z; or NULL. In the
  #  latter case the starting point is generated by method of moments.

  #--------------
  sub_Gomp_loglike_b = function(b, Z, theta1, theta2) {
    res = 0
    for (t in 2:T) {
      res = res + dnorm(
        Z[t], mean = -b * theta1 + (1 + b) * Z[t - 1],
        sd = sqrt(-b * (2 + b) * theta2), log = TRUE
      )
    }
    return(res)
  }
  #--------------


  ###############--------------- function script

  ### print start time if required
  if(!is.infinite(verbose)) {
    print(paste("GompPois_quasiGibbs: start time at ", Sys.time(), sep = ""))
  }

  ### check verbose
  if (verbose == 0) verbose = nsim + 1

  ### check burn
  if (burn <= 0) burn = 1




  ########## initialization

  ##### get starting values

  ### if starting values are not provided then use method of moments
  if (is.null(starter)) {
    starter = GompPois_MoM(Nstar)
    # check if estimated values belong to the support of the parameters
    if (starter$theta2 < 0.01) starter$theta2 = 0.01
    if (starter$b > -0.01) starter$b = -0.01
    if (starter$b < -0.99) starter$b = -0.99
    # sample Z from importance density
    starter$IS = c(GompPois_likelihood(
      Nstar, theta1 = starter$theta1, theta2 = starter$theta2, b = starter$b,
      nsim = 1e+4, log = TRUE, details = TRUE
    ))
    starter$Z = starter$IS$Z[, sample(
      c(1:1e+4), size = 1, replace = FALSE, prob = starter$IS$normWeights
    )]
  }

  ### initialize parameter values
  theta1 = starter$theta1
  theta2 = starter$theta2
  b = starter$b
  Z = starter$Z

  a = -b * theta1
  sigmaSq = -b * (2 + b) * theta2
  beta = log(-b / (1 + b))
  V = pi^2 / 3 # just for initialization, it is not used
  W = 1 # just for initialization, it is not used
  U = 1 # just for initialization, it is not used



  ### sample size
  T = length(Nstar)

  ### initialize inverse of B
  invB = matrix(0, nrow = T, ncol = T)

  ### storing vectors
  keep_theta1 = double(nsim)
  keep_theta2 = double(nsim)
  keep_b = double(nsim)




  ########## draw the chain
  for (insim in (1-burn):nsim) {

    ithin = 0

    ##### iterate up to thinning value
    while (ithin < thin) {

      ##########################################################################
      #                                                                        #
      #                           update U, V and W                            #
      #                                                                        #
      ##########################################################################

      U = rgamma(1, shape = phi1, rate = phi2 + 1 / theta2)
      V = 1 / rgamma(1, shape = 1, rate = 0.5 + beta^2 / (2 * c))
      W = 1 / rgamma(
        1, shape = 0.5 * (nu + 1),
        rate = 0.5 * (nu + (theta1 - zeta1)^2 / (theta2 * zeta2))
      )



      ##########################################################################
      #                                                                        #
      #                                update Z                                #
      #                                                                        #
      ##########################################################################

      for (t in 1:T) {

        if (t == 1) {
          mu = (theta1 * sigmaSq + (1 + b) * (Z[2] - a) * theta2) / (
            sigmaSq + (1 + b)^2 * theta2
          )
          tauSq = sigmaSq * theta2 / (sigmaSq + (1 + b)^2 * theta2)
        } else if (t == T) {
           mu = a + (1 + b) * Z[T-1]
           tauSq = sigmaSq
        } else {
           mu = (a + (1 + b) * (Z[t + 1] + Z[t-1] - a)) / (1 + (1 + b)^2)
           tauSq = sigmaSq / (1 + (1 + b)^2)
        }

        xi = Nstar[t] * tauSq + mu - lambertW_expArg(
          log(tauSq) + Nstar[t] * tauSq + mu
        )

        logC = -exp(xi) + xi * Nstar[t] - 0.5 / tauSq * (xi - mu)^2

        while (TRUE) {
          x = xi + sqrt(tauSq) * rnorm(100)
          check = -rexp(100) <= -exp(x) + x * Nstar[t] -
           0.5 / tauSq * (x - mu)^2 +
            0.5 / tauSq * (x - xi)^2 -
            logC
          if (sum(check) > 0) {
            Z[t] = x[check][1]
            break
          }
        }

      }



      ##########################################################################
      #                                                                        #
      #                                update b                                #
      #                                                                        #
      ##########################################################################

      strike = sub_Gomp_loglike_b(b, Z = Z, theta1 = theta1, theta2 = theta2) -
        rexp(1)
      delta = runif(1, min = 0, max = 2 * pi)
      delta_min = delta - 2 * pi
      delta_max = delta
      betaTilde = sqrt(c * V) * rnorm(1)

      while (TRUE) {
        betaProp = beta * cos(delta) + betaTilde * sin(delta)
        bProp = -1 / (1 + exp(-betaProp))
        loglike_bProp = sub_Gomp_loglike_b(
          bProp, Z = Z, theta1 = theta1, theta2 = theta2
        )
        if (loglike_bProp >= strike) {
          beta = betaProp
          b = bProp
          break
        } else {
          if (delta < 0) delta_min = delta else delta_max = delta
          delta = runif(1, min = delta_min, max = delta_max)
        }
      }



      ##########################################################################
      #                                                                        #
      #                        update theta1 and theta2                        #
      #                                                                        #
      ##########################################################################

      invB[1, 1] = 1
      invB[1, 2] = -(1 + b)
      invB[2, 1] = -(1 + b)
      invB[T, T] = 1
      for (t in 2:(T-1)) {
        invB[t, t] = 1 + (1 + b)^2
        invB[t, t + 1] = -(1 + b)
        invB[t + 1, t] = -(1 + b)
      }
      invB = invB / (1 - (1 + b)^2)

      commonFactor = 1 + W * zeta2 * sum(invB)
      rowSums_invB = c(rowSums(invB)) # it is B^-1 %*% 1_T
      #matrixQuadForm_theta2 = invB - W * zeta2 * tcrossprod(rowSums_invB) /
      #  commonFactor

      theta2 = 1 / rgamma(
        1, shape = 1 + 0.5 * T, rate = U + 0.5 * emulator::quad.form(
          invB - W * zeta2 * tcrossprod(rowSums_invB) / commonFactor, Z - zeta1
        )
      )

      theta1 = rnorm(
        1, mean = (W * zeta2 * crossprod(rowSums_invB, Z) + zeta1) /
          commonFactor,
        sd = sqrt(theta2 * W * zeta2 / commonFactor)
      )

      a = -b * theta1
      sigmaSq = -b * (2 + b) * theta2



      ##### end of single complete scan
      if (insim > 0) ithin = ithin + 1 else ithin = thin

    }

    ### check if outside of burn-in period
    if (insim > 0) {

      ### keep values
      keep_b[insim] = b
      keep_theta1[insim] = theta1
      keep_theta2[insim] = theta2

      ### print status of the chain
      if (insim %% verbose == 0) {
        print(paste(
          "iteration ", insim, " of ", nsim, " completed at time ", Sys.time(),
          sep = ""
        ))
      }

    }

  }

  ### print end time if required
  if (!is.infinite(verbose)) {
    print(paste("GompPois_quasiGibbs: end time at ", Sys.time(), sep = ""))
  }

  ### return results
  return(list(
    theta1 = keep_theta1,
    theta2 = keep_theta2,
    b = keep_b
  ))

}





############### appendix

if (FALSE) { # appendix script is not run

  # debugging
  {
    rm(list = ls()); gc(); cat("\14")
    setwd(paste(
      "/home/paolo/Desktop/PostDoc_Sapienza/VisitingToronto2023/",
      "AnimalMovement/codes/GMLossF-main_paolo/Functions/", sep = ""
    ))
    source("GompPois_ScaledPriors_quasiGibbs.R")
    source("extraFunctions.R")
    set.seed(1994)
    nsim = 1e+4
    T = 250
    theta1 = 1.9244
    theta2 = 0.4726 ^ 2 # that is 0.2233508
    b = -0.8 # -0.24

    Nstar = GompPois_simulate(T, theta1 = theta1, theta2 = theta2, b = b)
    rm(list = c("T", "theta1", "theta2", "b"))

    phi1 = 0.5
    phi2 = 0.5
    nu = 1
    zeta1 = 0
    zeta2 = 1
    c = 1
    starter = NULL
    burn = 1
    thin = 1
    verbose = 100
    start_time = Sys.time()
    res = GompPois_ScaledPriors_quasiGibbs(
      nsim, Nstar = Nstar, phi1 = phi1, phi2 = phi2, nu = nu, zeta1 = zeta1,
      zeta2 = zeta2, c = c, starter = starter, burn = burn, thin = thin,
      verbose = verbose
    )
    end_time = Sys.time()
  }

  difftime(end_time, start_time)

  mean(res$theta1)
  mean(res$theta2)
  mean(res$b)

  acf(res$theta1)
  acf(res$theta2)
  acf(res$b)
  #plot(density(res$b))

  plot(res$theta1, type = "l")
  plot(res$theta2, type = "l")
  plot(res$b, type = "l")

  sns::ess(res$theta1)
  sns::ess(res$theta2)
  sns::ess(res$b)

  ###############

  # check inversion of B

  b = -0.24
  T = 1e+3

  B = (1 + b)^(abs(outer(1 : T, 1 : T, "-"))) # get B matrix
  invB1 = chol2inv(chol(B))

  invB2 = matrix(0, nrow = T, ncol = T)
  invB2[1, 1] = 1
  invB2[1, 2] = -(1 + b)
  invB2[2, 1] = -(1 + b)
  invB2[T, T] = 1
  for (t in 2:(T-1)) {
    invB2[t, t] = 1 + (1 + b)^2
    invB2[t, t + 1] = -(1 + b)
    invB2[t + 1, t] = -(1 + b)
  }
  invB2 = invB2 / (1 - (1 + b)^2)

  invB1[1:3, 1:3]
  invB2[1:3, 1:3]

  B2 = chol2inv(chol(invB2))
  B2[1:3, 1:3]
  B[1:3, 1:3]

}